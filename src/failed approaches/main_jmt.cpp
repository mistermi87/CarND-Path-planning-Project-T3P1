#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Dense" //changed from core
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <typeinfo>
#include <cmath>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

//each trajectory is
struct Trajectory {
    vector<double> xt;
    vector<double> yt;
    double cost;
    double t_action_complete;
    };

int TIME =0;
vector<double> SIZE_OF_CARS_F {2, 1.2}; //length and width of car in frenet coordinates
double OVERALL_DT =0.02; //timestep between pathpoints
double A_MAX=3.5; //max. Acceleration m/s^2
double PREVIOUS_RECOG=5; //how many of the previous path's timepoints are taken into consideration for the new path
double MPH_IN_MS = 0.44704; //M/h in m/s
double VMAX=50*0.44704; //Vmax in m/s
double V_GOAL=45*0.44704; //V_goal in m/s
double S_FULL=6945.554;// length of the track
double CURRENT_LANE;
double LAST_LANE;

//Weights
double WEIGHT_DIST=1;   //How far does each trajectory take the car in 4s
double WEIGHT_SPEED=1;  //Does the car go over Vmax? is it going exceptionally slow?
double WEIGHT_ACCL=1;   //Dos the trajectory create high accelerations?
double WEIGHT_JERK=1;   //Dos the trajectory create jerk?
double WEIGHT_COLL=1;   //Does the trajectory lead to close distances to other cars or even collisions?


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<vector<double>> predictions_other_cars_xy(double time, vector<vector<double>> sensor_fusion){
    //[ id, x, y, vx, vy, s, d].
    vector<vector<double>> pred;
    for(int i=0;i<sensor_fusion.size();i++){
        vector<double> car_fut = {sensor_fusion[i][0],sensor_fusion[i][1]+time*sensor_fusion[i][3],sensor_fusion[i][2]+time*sensor_fusion[i][4]};
        pred.push_back(car_fut);
        }
    return pred;
}

vector<vector<double>> predictions_other_cars_sd(double time, vector<vector<double>> sensor_fusion, vector<double> map_waypoints_s, vector<double> map_waypoints_dx, vector<double> map_waypoints_dy ){
    //[ id, x, y, vx, vy, s, d].
    vector<vector<double>> pred;
    for(int i=0;i<sensor_fusion.size();i++){
        int prev_wp = -1;
        while(sensor_fusion[i][5] > map_waypoints_s[prev_wp+1] && (prev_wp < (int)(map_waypoints_s.size()-1) ))
        {
            prev_wp++;
            }
        int wp2 = (prev_wp+1)%map_waypoints_s.size();
        double d_fut=(map_waypoints_dx[prev_wp]*sensor_fusion[i][3]+map_waypoints_dy[prev_wp]*sensor_fusion[i][4])*time;
        double s_fut=(-map_waypoints_dy[prev_wp]*sensor_fusion[i][3]+map_waypoints_dx[prev_wp]*sensor_fusion[i][4])*time;
        vector<double> car_fut = {sensor_fusion[i][0],s_fut,d_fut};
        pred.push_back(car_fut);
        }
    return pred;
}

double closest_car_at_t(double time, double x, double y, vector<vector<double>> sensor_fusion){
    vector<vector<double>> other_cars_pos = predictions_other_cars_xy(time, sensor_fusion);
    double min_distance=9999;
    for(int i=0;i<other_cars_pos.size();i++){
        double dist_x=other_cars_pos[i][1]-x;
        double dist_y=other_cars_pos[i][2]-y;
        double dist_sq=dist_x*dist_x+dist_y*dist_y;
        if(dist_sq<min_distance){
            min_distance=dist_sq;
        }
    }
    return sqrt(min_distance);
}

vector<double> next_car_in_lane(vector<double> ego_car_data, vector<vector<double>> sensor_fusion){
    //[ id, x, y, vx, vy, s, d].
    double next_distance=99;
    int id=99;
    for(int i=0;i<sensor_fusion.size();i++){
        bool is_in_lane=abs(ego_car_data[3]-sensor_fusion[i][6])<2;
        double distance_to_car=fmod(sensor_fusion[i][5]-ego_car_data[2]+S_FULL,S_FULL);
        bool is_closest=distance_to_car<next_distance;
        bool is_in_front=distance_to_car>0;
        //cout<<sensor_fusion[i][0]<<" "<<distance_to_car<<" "<<abs(ego_car_data[3]-sensor_fusion[i][6])<< " is_in_lane "<< is_in_lane<<endl;
        if(is_in_lane & is_closest & is_in_front ){//is in same lane
                next_distance=distance_to_car;
                id=sensor_fusion[i][0];

        }
    }
    //cout<<"closest: "<<id<<" Distance: "<<next_distance<<endl;
    double speed=0;
    if(id<99){
        speed=sqrt(sensor_fusion[id][3]*sensor_fusion[id][3]+sensor_fusion[id][4]*sensor_fusion[id][4]);
    }
    return {id,next_distance,speed};
}


vector<double> JMT(vector< double> start, vector <double> end, double T)
{
    double s_i=start[0];
    double s_i_dot=start[1];
    double s_i_double_dot=start[2];

    double s_f=end[0];
    double s_f_dot=end[1];
    double s_f_double_dot=end[2];

    double a_0=s_i;
    double a_1=s_i_dot;
    double a_2=0.5*s_i_double_dot;


    MatrixXd T_matrix = MatrixXd(3, 3);
    double T2=T*T;
    double T3=T*T*T;
    double T4=T*T*T*T;
    double T5=T*T*T*T*T;

    T_matrix<<  T3,     T4,     T5,
                3*T2,   4*T3,   5*T4,
                6*T,   12*T2,  20*T3;

    VectorXd solved_side = VectorXd(3);
    solved_side << s_f - (s_i + s_i_dot*T + 0.5*s_i_double_dot*T*T),
                   s_f_dot - (s_i_dot + s_i_double_dot*T),
                   s_f_double_dot - s_i_double_dot;

    VectorXd a_345;

    MatrixXd T_matrix_inv=T_matrix.inverse();

    a_345=T_matrix_inv*solved_side;

    double a_3=a_345[0];
    double a_4=a_345[1];
    double a_5=a_345[2];


    vector<double> result = {a_0, a_1, a_2, a_3, a_4, a_5};

//    cout<<endl;
//	for(int i = 0; i < 6; i++)
//	{
//	    cout<<result[i]<<" ";
//	}
//	cout<<endl;

    return result;
}

double get_p_value(vector<double> poly, double x){
    double mult=1;
    double result=0;
    for(int i=0; i<poly.size();i++){
        result+=poly[i]*mult;
        mult*=x;
    }
    return result;
}

double get_qp_grad1(vector<double> poly, double x){
    double mult=1;
    double result=0;
    for(int i=0; i<poly.size()-1;i++){
        result+=poly[i+1]*mult*(i+1);
        mult*=x;
    }
    return result;
}

double get_qp_grad2(vector<double> poly, double x){
    double mult=1;
    double result=0;
    for(int i=0; i<poly.size()-2;i++){
        result+=poly[i+2]*mult*(i+2)*(i+1);
        mult*=x;
    }
    return result;
}

double get_qp_grad3(vector<double> poly, double x){
    double mult=1;
    double result=0;
    for(int i=0; i<poly.size()-3;i++){
        result+=poly[i+3]*mult*(i+3)*(i+2)*(i+1);
        mult*=x;
    }
    return result;
}


Trajectory create_splines_xt_yt(vector<double> car_data,vector<double> previous_x, vector<double> previous_y, vector<double> goal, vector<double> map_waypoints_s, vector<double> map_waypoints_x, vector<double> map_waypoints_y, vector<vector<double>> sensor_fusion){


    //Previous path+10
    int path_size=previous_x.size();
    double t_p=path_size*OVERALL_DT;
    double x_p=previous_x[path_size-1];
    double y_p=previous_y[path_size-1];

    double x_p_m1=previous_x[path_size-2];
    double y_p_m1=previous_y[path_size-2];

    double x_p_m2=previous_x[path_size-3];
    double y_p_m2=previous_y[path_size-3];

    double x_p_m3=previous_x[path_size-4];
    double y_p_m3=previous_y[path_size-4];

    double v_x_p=(x_p-x_p_m1)/OVERALL_DT;
    double v_y_p=(y_p-y_p_m1)/OVERALL_DT;

    double v_x_p_m1=(x_p_m1-x_p_m2)/OVERALL_DT;
    double v_y_p_m1=(y_p_m1-y_p_m2)/OVERALL_DT;

    double v_x_p_m2=(x_p_m2-x_p_m3)/OVERALL_DT;
    double v_y_p_m2=(y_p_m2-y_p_m3)/OVERALL_DT;

    double a_x_p=(v_x_p-v_x_p_m1)/OVERALL_DT;
    double a_y_p=(v_y_p-v_y_p_m1)/OVERALL_DT;

    double a_x_p_m1=(v_x_p_m1-v_x_p_m2)/OVERALL_DT;
    double a_y_p_m1=(v_y_p_m1-v_y_p_m2)/OVERALL_DT;

    cout<<"0: X: " <<x_p<<" Y: "<<y_p<<" X-1: "<<x_p_m1<<" Y-1: "<<y_p_m1<<" X-2: "<<x_p_m2<<" Y-2: "<<y_p_m2<<" Vx: "<<v_x_p<<" Vy: "<<v_y_p<<" Vx-1: "<<v_x_p_m1<<" Vy-1: "<<v_y_p_m1<< endl;

    //--------------------------------------------------------------------------------

    // at goal
    //where you want to get
    // assumptions: a_d=0, v_d0=0;



    vector<double> goalXY= getXY(goal[0], goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);
    vector<double> goalXY2= getXY(goal[0]+2, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);//second point for getting direction at goal (along s)

    double gradient_x_goal=(goalXY2[0]-goalXY[0])/distance(goalXY[0],goalXY[1],goalXY2[0],goalXY2[1]);
    double gradient_y_goal=(goalXY2[1]-goalXY[1])/distance(goalXY[0],goalXY[1],goalXY2[0],goalXY2[1]);

    double goal_velocity=goal[2];
    double current_velocity=car_data[5];

    //double a_d_factor=(goal_velocity-current_velocity)/abs(goal_velocity-current_velocity);

    double distance_to_goal=distance(x_p,y_p,goalXY[0],goalXY[1]); //approximation of length of spline from current position to goal (the real result would be an Integration of distances of x(t),y(t)) maybe estimated with "Gaussian quadrature")
    double time_to_goal_vel = abs((goal_velocity-current_velocity)/A_MAX);
    double distance_to_goal_vel=(0.5*current_velocity+0.5*goal_velocity)*time_to_goal_vel;
    double t_1=99;
    double a_s_1=0;
    if(distance_to_goal<distance_to_goal_vel){
        a_s_1=(goal_velocity-current_velocity)/abs(goal_velocity-current_velocity)*A_MAX; //positive acceleration if goal velocity is not jet reached
        goal_velocity=sqrt(current_velocity*current_velocity+2*distance_to_goal*a_s_1);
        t_1=abs((goal_velocity-current_velocity)/A_MAX);

    }
    else{
        t_1=time_to_goal_vel+(distance_to_goal-distance_to_goal_vel)/goal_velocity;
    }
    //cout<<"goal_velocity "<<goal_velocity<<" current_velocity: "<<current_velocity<<" distance_to_goal "<<distance_to_goal<<" time_to_goal_vel "<< time_to_goal_vel <<" distance_to_goal_vel "<<distance_to_goal_vel <<" t_1 "<<t_1<<endl;
    double v_x_1=gradient_x_goal*goal_velocity;
    double v_y_1=gradient_y_goal*goal_velocity;

    double a_x_1=gradient_x_goal*a_s_1;
    double a_y_1=gradient_y_goal*a_s_1;

    double x_1=goalXY[0];
    double y_1=goalXY[1];

    //cout<<"1:    "<<t_1<<" X1: " <<x_1<<" Y1:  "<<y_1<<" Vx1: "<<v_x_1<<" Vy1: "<<v_y_1<<" Ax1: "<<a_x_1<<" Ay1: "<<a_y_1<< endl;



    //at horizon
    //end of spline. All trajectories are judged (via cost function) up to this point (~80m away at opt. speed) as most actions can be completed within this distance

    double time_horizon=4; //in sec
    //at low velocities the goal might be further away than 4 sec
    if(time_horizon<t_1){
        time_horizon=t_1+OVERALL_DT;
    }
    double t_accl=0;
    if(goal[2]!=goal_velocity){
        t_accl=abs((goal[2]-goal_velocity)/A_MAX);
    }
    double distance_horizon=goal_velocity*t_accl+ 0.5*A_MAX*t_accl*t_accl+ goal[2]*((time_horizon-t_1)-t_accl);

    vector<double> horizon= getXY(goal[0]+distance_horizon, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);
    vector<double> horizon2= getXY(goal[0]+distance_horizon+1, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);

    double gradient_x_horizon=(horizon2[0]-horizon[0])/distance(horizon[0],horizon[1],horizon2[0],horizon2[1]);
    double gradient_y_horizon=(horizon2[1]-horizon[1])/distance(horizon[0],horizon[1],horizon2[0],horizon2[1]);

    double t_2=time_horizon;
    double x_2=horizon[0];
    double y_2=horizon[1];

    double v_x_2=gradient_x_horizon*goal[2];
    double v_y_2=gradient_x_horizon*goal[2];

    //cout<<"Hori: "<<t_2<<" " <<x_2<<" "<<y_2<<endl;


    //polynome from start to goal
    vector<double> p_1_xt=JMT({x_p,(v_x_p+v_x_p_m1)/2,(a_x_p+a_x_p_m1)/2}, {x_1,v_x_1,a_x_1}, t_1);
    vector<double> p_1_yt=JMT({y_p,(v_y_p+v_y_p_m1)/2,(a_y_p+a_y_p_m1)/2}, {y_1,v_y_1,a_y_1}, t_1);

    //polynome from goal to horizon
    vector<double> p_2_xt=JMT({x_1,v_x_1,a_x_1}, {x_2,v_x_2,0}, time_horizon-t_1);
    vector<double> p_2_yt=JMT({y_1,v_y_1,a_y_1}, {y_2,v_y_2,0}, time_horizon-t_1);


    //calculate cost of trajectory
    //==========================================================================
    double cost_speed=0;
    double cost_accl=0;
    double cost_jerk=0;
    double cost_collision=0;

    double last_acc_x =0;
    double last_acc_y =0;

    double time_step=OVERALL_DT*5;

    for(double time=0; time>(t_1-t_p);time+=time_step){
        // speed above Vmax?
        double x_t=get_p_value(p_1_xt,time);
        double y_t=get_p_value(p_1_yt,time);
        double v_x_t=get_qp_grad1(p_1_xt,time);
        double v_y_t=get_qp_grad1(p_1_yt,time);
        double a_x_t=get_qp_grad2(p_1_xt,time);
        double a_y_t=get_qp_grad2(p_1_yt,time);
        double j_x_t=get_qp_grad3(p_1_xt,time);
        double j_y_t=get_qp_grad3(p_1_yt,time);

        double speed_rel_to_vmax=VMAX-sqrt(v_x_t*v_x_t+v_y_t+v_y_t);
        if(speed_rel_to_vmax<0) {speed_rel_to_vmax = 0;} //only count positive values (fast) https://stackoverflow.com/questions/11275187/r-replacing-negative-values-by-zero
        cost_speed+=speed_rel_to_vmax;

        //acceleration
        cost_accl+=a_x_t*a_x_t+a_y_t*a_y_t;

        // Jerk?
        cost_jerk+=j_x_t*j_x_t+j_y_t*j_y_t;

        // Collision


        double closest_car=4-closest_car_at_t(time+t_p, x_t, y_t, sensor_fusion);
        double c3=-closest_car*closest_car*closest_car;
        if(c3<0) {c3 = 0;} //only count positive values (fast) (see above)
        cost_collision+=c3;
        }

    // will there be slow cars in front of vehicle on possible new lane?
    // will there be jerk or acceleration necessary to stay on final lane
    time_step=OVERALL_DT*20;

    for(double time=0; time>(time_horizon-t_1);time+=time_step){
        // speed above Vmax?
        double x_t=get_p_value(p_2_xt,time);
        double y_t=get_p_value(p_2_yt,time);
        double v_x_t=get_qp_grad1(p_2_xt,time);
        double v_y_t=get_qp_grad1(p_2_yt,time);
        double a_x_t=get_qp_grad2(p_2_xt,time);
        double a_y_t=get_qp_grad2(p_2_yt,time);
        double j_x_t=get_qp_grad3(p_2_xt,time);
        double j_y_t=get_qp_grad3(p_2_yt,time);

        double speed_rel_to_vmax=VMAX-sqrt(v_x_t*v_x_t+v_y_t+v_y_t);
        if(speed_rel_to_vmax<0) {speed_rel_to_vmax = 0;} //only count positive values (fast) https://stackoverflow.com/questions/11275187/r-replacing-negative-values-by-zero
        cost_speed+=speed_rel_to_vmax;

        //acceleration
        cost_accl+=a_x_t*a_x_t+a_y_t*a_y_t;

        // Jerk?
        cost_jerk+=j_x_t*j_x_t+j_y_t*j_y_t;

        // Collision only until goal
        double closest_car=4-closest_car_at_t(time+t_1, x_t, y_t, sensor_fusion);
        double c3=closest_car*closest_car*closest_car;
        if(c3<0) {c3 = 0;} //only count positive values (fast) (see above)
        cost_collision+=c3;
        }

    // how far in 4s (s)
    double cost_distance=exp(-0.1*(goal[0]+distance_horizon));


    Trajectory this_traj;
    this_traj.xt=p_1_xt;
    this_traj.yt=p_1_yt;
    this_traj.cost=WEIGHT_SPEED*cost_speed+WEIGHT_ACCL*cost_accl+WEIGHT_JERK*cost_jerk+WEIGHT_COLL*cost_collision+WEIGHT_DIST*cost_distance;
    this_traj.t_action_complete=t_1;
    return this_traj;
}





int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];



          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	//int prev_path_size=previous_path_x.size();


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds


            vector<double> ego_car_data{car_x,car_y,car_s,car_d,car_yaw,car_speed*0.44704}; //<-- speed conversion to m/s;
            cout<<"--- --- ---"<<endl;
            cout<<" X: "<<ego_car_data[0]<<" Y: "<<ego_car_data[1]<<" S: "<<ego_car_data[2]<<" D: "<<ego_car_data[3]<<" Yaw: "<<ego_car_data[4]<<" v: "<<ego_car_data[5]<<endl;

          	 //get info on road for spline of road around ego-vehicle using waypoints where 0 is closest waypoint

            //int next_wp = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
            //int aft_next_wp=next_wp%map_waypoints_x.size();




          	//check distances to surrounding cars
          	//if lane is clear --> constant speed
          	vector<double> goal;
          	vector<double> ncil=next_car_in_lane(ego_car_data, sensor_fusion);
          	cout<<ncil[0]<<" "<<ncil[1]<<" "<<ncil[2]<<endl;

            if((ncil[0]<99) & (ncil[1]<40)){
                goal={car_s+ncil[1]-10, 6, ncil[2]};
                cout<<"next car is close"<<endl;
            }else{
                goal={car_s+50, 6, 20};
                cout<<"next car is far"<<endl;
            }



          	//if car is in front make decision and rate cost


          double pos_x;
          double pos_y;
          double angle;
          int path_size = previous_path_x.size();

          if(path_size<(PREVIOUS_RECOG)){
            for(int i=0; i<PREVIOUS_RECOG-path_size;i++){
                previous_path_x.insert(previous_path_x.begin(),car_x);
                previous_path_y.insert(previous_path_y.begin(),car_y);
            }
          }

          for(int i = 0; i <previous_path_x.size() ; i++)
          {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
          }
//          cout<<endl;
//
//          for(int i=0; i<path_size;i++){
//                cout<<previous_path_x[i]<<" \t "<<previous_path_y[i]<<endl;
//                }
//          cout<<endl;

          //vector<double> goal{car_s+50, 6, 20};
          Trajectory next_traj=create_splines_xt_yt(ego_car_data, previous_path_x, previous_path_y, goal,  map_waypoints_s,  map_waypoints_x, map_waypoints_y, sensor_fusion );

          //cout<<"Next Step:" <<  xtyt[0](OVERALL_DT) <<" " << xtyt[1](OVERALL_DT)<<endl;


          double how_much_to_add= std::fmin(abs(50-previous_path_x.size()),(next_traj.t_action_complete/OVERALL_DT)); //maximum of 50 points to the path
          cout<<how_much_to_add<< " new points added"<<endl;

          for(int i = 0; i < how_much_to_add; i++)
          {
              next_x_vals.push_back(get_p_value(next_traj.xt,OVERALL_DT*i));
              next_y_vals.push_back(get_p_value(next_traj.yt,OVERALL_DT*i));
          }

          //END

          json msgJson;

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
