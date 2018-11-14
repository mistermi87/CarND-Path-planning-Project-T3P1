#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
//#include "Functions.cpp"
#include "spline.h"
#include <typeinfo>

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


struct Trajectory {
    tk::spline xt;
    tk::spline yt;
    double cost;
    double t_action_complete;
    bool includes_lanechange;
    };

struct Cars_around {
    int id_infront;
    int id_r;
    int id_l;
    int id_r_infront;
    int id_l_infront;
    double distance;
    double distance_r;
    double distance_l;
    double speed;
    double speed_r;
    double speed_l;
    };

double TIME =0;
double OVERALL_DT =0.02; //timestep between pathpoints
double A_MAX=3; //max. Acceleration m/s^2
double MPH_IN_MS = 0.44704; //M/h in m/s
double VMAX=47*MPH_IN_MS; //Vmax in m/s
//double V_GOAL=VMAX; //V_goal in m/s
double S_FULL=6945.554;// length of the track
double CURRENT_LANE=1;
double GOAL_LANE=1;

double GOAL_TIME=0; //need this?
bool IS_CHANGING_LANES;
double KEEP_DISTANCE=20; //~the distance around which the ego-car regulates its speed to match up with car in front


//Weights
double WEIGHT_DIST=100;   //How far does each trajectory take the car in 4s
double WEIGHT_SPEED=1;  //Does the car go over Vmax? is it going exceptionally slow?
double WEIGHT_ACCL=100;   //Dos the trajectory create high accelerations?
double WEIGHT_JERK=100;   //Dos the trajectory create jerk?
double WEIGHT_COLL=10000;   //Does the trajectory lead to close distances to other cars or even collisions?


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
        vector<double> car_fut={sensor_fusion[i][0],sensor_fusion[i][1]+time*sensor_fusion[i][3],sensor_fusion[i][2]+time*sensor_fusion[i][4]};
        pred.push_back(car_fut);
        }
    return pred;
}

//vector<vector<double>> predictions_other_cars_sd(double time, vector<vector<double>> sensor_fusion, vector<double> map_waypoints_s, vector<double> map_waypoints_dx, vector<double> map_waypoints_dy ){
//    //[ id, x, y, vx, vy, s, d].
//    vector<vector<double>> pred;
//    vector<double> car_fut;
//    for(int i=0;i<sensor_fusion.size();i++){
//        int prev_wp = -1;
//        while(sensor_fusion[i][5] > map_waypoints_s[prev_wp+1] && (prev_wp < (int)(map_waypoints_s.size()-1) ))
//        {
//            prev_wp++;
//        }
//        int wp2 = (prev_wp+1)%map_waypoints_s.size();
//        double d_fut=(map_waypoints_dx[prev_wp]*sensor_fusion[i][3]+map_waypoints_dy[prev_wp]*sensor_fusion[i][4])*time;
//        double s_fut=(-map_waypoints_dy[prev_wp]*sensor_fusion[i][3]+map_waypoints_dx[prev_wp]*sensor_fusion[i][4])*time;
//        car_fut={sensor_fusion[i][0],s_fut,d_fut};
//        pred.push_back(car_fut);
//        }
//    return pred;
//}

double closest_car_at_t(double time, double x, double y, vector<vector<double>> sensor_fusion){
    vector<vector<double>> other_cars_pos=predictions_other_cars_xy(time, sensor_fusion);
    double min_distance=9999;
    //double id;
    for(int i=0;i<other_cars_pos.size();i++){
        double dist_x=other_cars_pos[i][1]-x;
        double dist_y=other_cars_pos[i][2]-y;
        double dist_sq=dist_x*dist_x+dist_y*dist_y;
        if(dist_sq<min_distance){
            min_distance=dist_sq;
            //id=sensor_fusion[i][0];
        }
    }
    return sqrt(min_distance);
}

Cars_around next_car_closeby(vector<double> ego_car_data, vector<vector<double>> sensor_fusion){
    //[ id, x, y, vx, vy, s, d].
    //computes the next car in the same lane, as well as whether there is a car directly left or right of the ego car
    Cars_around cars;
    cars.distance=999;
    cars.distance_l=999;
    cars.distance_r=999;
    cars.id_infront=99;
    cars.id_l=99;
    cars.id_r=99;
    cars.id_l_infront=99;
    cars.id_r_infront=99;
    cars.speed=0;
    cars.speed_l=0;
    cars.speed_r=0;

    double my_d=ego_car_data[3];
    double my_s=ego_car_data[2];
    for(int i=0;i<sensor_fusion.size();i++){
        double car_d=sensor_fusion[i][6];
        double car_s=sensor_fusion[i][5];


        double distance_to_car=car_s-my_s;//fmod(car_s-my_s+S_FULL,S_FULL);//
        bool is_in_front=distance_to_car>3;


        bool is_in_lane=abs(my_d-car_d)<2;
        bool is_closest_in_lane=(is_in_lane & distance_to_car<cars.distance);

        //identify car in front of me
        if(is_closest_in_lane & is_in_front ){
                cars.distance=distance_to_car;
                cars.id_infront=sensor_fusion[i][0];
        }


        bool is_parallel_to_me=((distance_to_car<10) & (distance_to_car>-12));



        //identifying car left of me
        bool is_left= ((3< (my_d-car_d)) &  ((my_d-car_d)<6));
        if(is_parallel_to_me & is_left){
            cars.id_l=sensor_fusion[i][0];
        }
        //identifying car right of me
        bool is_right= ((-3> (my_d-car_d)) &  ((my_d-car_d) >-6));
        if(is_parallel_to_me & is_right){
            cars.id_r=sensor_fusion[i][0];
        }

        //identifying car in front of me in left lane
        bool is_closest_on_left=(is_left & distance_to_car<cars.distance_l);
        if(is_closest_on_left & is_in_front ){//is in same lane
                cars.distance_l=distance_to_car;
                cars.id_l_infront=sensor_fusion[i][0];
        }

        //identifying car in front of me in right lane
        bool is_closest_on_right=(is_right & distance_to_car<cars.distance_r);
        if(is_closest_on_right & is_in_front ){//is in same lane
                cars.distance_r=distance_to_car;
                cars.id_r_infront=sensor_fusion[i][0];
        }


    }
    //cout<<"closest: "<<id<<" Distance: "<<next_distance<<endl;
    if(cars.id_infront<99){
        cars.speed=sqrt(sensor_fusion[cars.id_infront][3]*sensor_fusion[cars.id_infront][3]+sensor_fusion[cars.id_infront][4]*sensor_fusion[cars.id_infront][4]);
    }
    if(cars.id_l_infront<99){
        cars.speed_l=sqrt(sensor_fusion[cars.id_l_infront][3]*sensor_fusion[cars.id_l_infront][3]+sensor_fusion[cars.id_l_infront][4]*sensor_fusion[cars.id_l_infront][4]);
    }
    if(cars.id_r_infront<99){
        cars.speed_r=sqrt(sensor_fusion[cars.id_r_infront][3]*sensor_fusion[cars.id_r_infront][3]+sensor_fusion[cars.id_r_infront][4]*sensor_fusion[cars.id_r_infront][4]);
    }

    //cout<<"l"<<'\t'<<"c"<<'\t'<<"r"<<endl;
    //cout<<cars.id_l_infront<<'\t'<<cars.id_infront<<'\t'<<cars.id_r_infront<<endl;
    //cout<<cars.speed_l<<'\t'<<"x"<<'\t'<<cars.speed_r<<endl;
    //cout<<cars.id_l<<'\t'<<"x"<<'\t'<<cars.id_r<<endl;

    return cars;
}

Trajectory create_splines_xt_yt(vector<double> car_data,vector<double> previous_x, vector<double> previous_y, vector<double> goal, vector<double> map_waypoints_s, vector<double> map_waypoints_x, vector<double> map_waypoints_y, vector<vector<double>> sensor_fusion){
    //Creates two splines x(t) and y(t) that make up the trajectory. This trajectory's  quality is then measured with several cost function

    int path_size=previous_x.size();


    /*
    STARTING POINTS
    */

    double t_m2;
    double x_m2;
    double y_m2;
    double t_m1;
    double x_m1;
    double y_m1;
    double t_0;
    double x_0;
    double y_0;


    if(path_size>2){

        // Two before Previous path end
        t_m2=OVERALL_DT*(path_size-2);
        x_m2=previous_x[path_size-3];
        y_m2=previous_y[path_size-3];

        // One before Previous path end
        t_m1=OVERALL_DT*(path_size-1);
        x_m1=previous_x[path_size-2];
        y_m1=previous_y[path_size-2];

        //Previous path end
        t_0=OVERALL_DT*path_size;
        x_0=previous_x[path_size-1];
        y_0=previous_y[path_size-1];

        //--> three points are used in order approximate current acceleration in spline


    }else{
        //for the first (few) frames
        t_m2=-2;
        x_m2=car_data[0] - 2* OVERALL_DT * cos(car_data[4]);
        y_m2=car_data[1] - 2* OVERALL_DT * sin(car_data[4]);


        t_m1=-1;
        x_m1=car_data[0] - OVERALL_DT *cos(car_data[4]);
        y_m1=car_data[1] - OVERALL_DT *sin(car_data[4]);

        t_0=0;
        x_0=car_data[0];
        y_0=car_data[1];



    }

    /*
    GOAL POINTS
    */
    //where car should end up
    // a second point is introduced to force let the splines towards a certain gradient (goal velocity)
    // an accelarion of 0 is assumed at the goal (simplification)


    //Two-Step approximation of length of spline from current position to goal (the real result would be an Integration of distances of x(t),y(t)) maybe estimated with "Gaussian quadrature")

    vector<double> halfway_point;
    if(car_data[2]<goal[0]){
        halfway_point=getXY((goal[0]+car_data[2])/2, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);
    }else{
        halfway_point=getXY(fmod((goal[0]+S_FULL+car_data[2])/2,S_FULL), goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);//after one round the goal is in the second lap (eg s=5) while the car is still in the first lap (s=~6800)
    }
    vector<double> goalXY= getXY(goal[0], goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);
    double distance_to_goal=distance(x_0,y_0,halfway_point[0],halfway_point[1])+distance(halfway_point[0],halfway_point[1],goalXY[0],goalXY[1]);



    double goal_velocity=goal[2];
    double current_velocity=car_data[5];

    double time_to_goal_vel = abs(goal_velocity-current_velocity)/A_MAX;
    double distance_to_goal_vel=(0.5*current_velocity+0.5*goal_velocity)*time_to_goal_vel;
    double time_between_goal_velocity_and_goal=(distance_to_goal-distance_to_goal_vel)/(0.5*current_velocity+0.5*goal_velocity);// /goal_velocity;


    double t_1;
    if(distance_to_goal_vel>distance_to_goal){
        goal[0]=goal[0]-distance_to_goal+distance_to_goal_vel;
        goalXY= getXY(goal[0], goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);
        t_1=time_to_goal_vel+t_0; //if goal velocity can not be reached until the goal distance this term comes into action. It is a rough approximation and could introduce a too high acceleration
        cout<<goal[0]<<" "<<t_1-t_0<<endl;
    }else{
        t_1=time_to_goal_vel+time_between_goal_velocity_and_goal+t_0;
    }

    double x_1=goalXY[0];
    double y_1=goalXY[1];

    //cout<<"goal_velocity "<<goal_velocity<<" current_velocity: "<<current_velocity<<" distance_to_goal "<<distance_to_goal<<" t_1 "<<t_1<<endl;


    //second point after goal to set the yaw and the velocity around the goal point
    double goal_2_s=goal[0]+OVERALL_DT*4*goal_velocity;
    vector<double> goal_2_XY= getXY(goal_2_s, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);

    double t_2=t_1+OVERALL_DT*4;
    double x_2=goal_2_XY[0];
    double y_2=goal_2_XY[1];

    //cout<<"2:    "<<t_2<<" " <<x_2<<" "<<y_2<<endl;

    /*
    HORIZON
    */
    //end of spline. All trajectories distance is judged (via cost function) up to this point ( 5 Seconds or ~100m away) as most actions can be completed within this distance.

    double time_horizon=5; //in sec

    //at low velocities the goal might be further away than that time
    if(time_horizon<t_2){
        time_horizon=t_2+OVERALL_DT*4;
    }

    double distance_horizon=goal_velocity*(time_horizon-t_1);

    vector<double> horizon= getXY(goal[0]+distance_horizon, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y);

    double t_3=time_horizon;
    double x_3=horizon[0];
    double y_3=horizon[1];



    tk::spline xt;
    xt.set_points({t_m2,t_m1,t_0,t_1,t_2,t_3},{x_m2,x_m1,x_0,x_1,x_2,x_3});


    tk::spline yt;
    yt.set_points({t_m2,t_m1,t_0,t_1,t_2,t_3},{y_m2,y_m1,y_0,y_1,y_2,y_3});

    //calculate cost of trajectory
    //==========================================================================
    double cost_speed=0;
    double cost_accl=0;
    double cost_jerk=0;
    double cost_collision=0;

    //stepsize for calculating cost
    double time_step=OVERALL_DT*5;

    //starting with approximation values
    double last_x;
    double last_y;
    double last_v_x;
    double last_v_y;
    double last_acc_x;
    double last_acc_y;
    cout<<"here"<<endl;

    if(path_size>6){ //step-size-1
        last_x=previous_x[path_size-6];
        last_y=previous_y[path_size-6];
        last_v_x=(previous_x[path_size-5]-last_x)/OVERALL_DT;
        last_v_y=(previous_y[path_size-5]-last_y)/OVERALL_DT;
        last_acc_x=((previous_x[path_size-1]-previous_x[path_size-2])/OVERALL_DT-last_v_x)/(4*OVERALL_DT);
        last_acc_y=((previous_y[path_size-1]-previous_y[path_size-2])/OVERALL_DT-last_v_y)/(4*OVERALL_DT);

    }else{
        last_x=car_data[0];
        last_y=car_data[1];
        last_v_x=0;
        last_v_y=0;
        last_acc_x=0;
        last_acc_y=0;
    }


//iterating over data points
    for(float time=t_0; time<t_1;time+=time_step){
        // speed above Vmax?
        double x=xt(time);
        double y=yt(time);
        double v_x=(x-last_x)/time_step;
        double v_y=(y-last_y)/time_step;
        double acc_x=(v_x-last_v_x)/time_step;
        double acc_y=(v_y-last_v_y)/time_step;
        double jerk_x=(acc_x-last_acc_x)/time_step;
        double jerk_y=(acc_y-last_acc_y)/time_step;
        //cout<<"V "<<v_x<<" "<<v_y<<" A "<<acc_x<<" "<<acc_y<<" || ";

        //speed
        double speed_to_vmax=VMAX-sqrt(v_x*v_x+v_y*v_y);
        if(speed_to_vmax<0) {speed_to_vmax = 0;} //only count positive values
        cost_speed+=speed_to_vmax;

        //acceleration
        double accl_tot=acc_x*acc_x+acc_y*acc_y; //omitted sqrt, unnecessary computation, can be taken into account in the following steps
        if(accl_tot<100) {accl_tot = 0;}
        cost_accl+=accl_tot;

        // Jerk?
        double jerk_tot=jerk_x*jerk_x+jerk_y*jerk_y; //omitted sqrt, unnecessary computation, can be taken into account in the following steps
        if(jerk_tot<100) {jerk_tot = 0;}
        cost_jerk+=jerk_tot;


        //Collision other car, comes within 1m to our trajectory (low value is neeeded otherwise cars of other lanes will be extrapolated (using v_x,v_y=const. in "closest_car_at_t") into collision course when cars are in a curve
        double closest_car=closest_car_at_t(time, x, y, sensor_fusion);
        //double c3=closest_car*closest_car*closest_car;
        if(closest_car<1) {
            closest_car = 1;
        } else {
            closest_car = 0;
            }
        cost_collision+=exp(-time)*closest_car ;//weight by distance along path

        last_x=x;
        last_y=y;
        last_v_x=v_x;
        last_v_y=v_y;
        last_acc_x =acc_x;
        last_acc_y =acc_y;
        }

    // how far in 4s (s)
    double cost_distance=1000*exp(-0.1*(goal[0]+distance_horizon-car_data[2]));

    cout<<"Speed: "<<cost_speed<<" Acc: "<<cost_accl<<" Jerk: "<<cost_jerk<<" Dist: "<<cost_distance<<" Collision: "<<cost_collision<<endl;
    Trajectory this_traj;
    this_traj.xt=xt;
    this_traj.yt=yt;
    this_traj.cost=WEIGHT_COLL*cost_collision+WEIGHT_DIST*cost_distance+WEIGHT_SPEED*cost_speed+WEIGHT_JERK*cost_jerk+WEIGHT_ACCL*cost_accl;
    this_traj.t_action_complete=t_1;
    this_traj.includes_lanechange=abs(goal[1]-car_data[3])>2;
    //cout<<"lane-change distance"<<goal[1]-car_data[3]<<endl;
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

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	int path_size = previous_path_x.size();

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds


            vector<double> ego_car_data{car_x,car_y,car_s,car_d,car_yaw,car_speed*0.44704}; //<-- speed conversion to m/s;
            //cout<<"--- --- ---"<<endl;
            cout<<" X: "<<ego_car_data[0]<<" Y: "<<ego_car_data[1]<<" S: "<<ego_car_data[2]<<" D: "<<ego_car_data[3]<<" Yaw: "<<ego_car_data[4]<<" v: "<<ego_car_data[5]<<endl;

          	 //get info on road for spline of road around ego-vehicle using waypoints where 0 is closest waypoint

            //int next_wp = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
            //int aft_next_wp=next_wp%map_waypoints_x.size();


            CURRENT_LANE=round((ego_car_data[3]-2)/4);

            vector<vector<double>> goals;
            double far_goal = 80;
            double near_goal = 50;


            //check distances to surrounding cars and categorize them for decision making
          	Cars_around cars_closeby=next_car_closeby(ego_car_data, sensor_fusion);


          	//check distances to surrounding cars for emergency braking
//            bool crash_now = closest_car_at_t(0,ego_car_data[0],ego_car_data[1],sensor_fusion)<3;
//
//            //emergency braking
//            if(crash_now){
//
//                //double min_distance=path_size*OVERALL_DT*VMAX;
//                //double min_distance_braking=0.5*VMAX*VMAX/A_MAX;
//                goals.push_back({car_s+10, CURRENT_LANE*4+2, 1});
//                cout<<"EMERGENCY BRAKE"<<endl;
//
//            }else{

                goals.push_back({car_s+near_goal, GOAL_LANE*4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance/5))}); //if no car in front --> distance=99 --> goal v= ~Vmax;

                //Only decide while in lane
                if(IS_CHANGING_LANES==0){// | crash_now | crash_soon | cars_closeby.distance<3 | cars_closeby.distance_l<3 | cars_closeby.distance_r<3){
                    //If car in front is detected, alternatives are calculated
                    if((cars_closeby.id_infront<99) & (cars_closeby.distance<40)){
                        //when in right lane AND there is no car left of me AND ( there is no car on front of me on the left lane OR the car in the left lane will be further in front in 4 seconds than the car in this lane)
                        if(CURRENT_LANE==2 & cars_closeby.id_l==99 & (cars_closeby.id_l_infront==99 | (cars_closeby.distance_l+cars_closeby.speed_l*4)>(cars_closeby.distance+cars_closeby.speed*4))){
                            goals.push_back({car_s+near_goal, 4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_l/5))});
                            goals.push_back({car_s+far_goal, 4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_l/5))});
                            }
                        //when in center lane...
                        if(CURRENT_LANE==1){
                            if(cars_closeby.id_l==99 & (cars_closeby.id_l_infront==99 | (cars_closeby.distance_l+cars_closeby.speed_l*4)>(cars_closeby.distance+cars_closeby.speed*4))){
                                goals.push_back({car_s+near_goal, 2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_l/5))});
                                goals.push_back({car_s+far_goal, 2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_l/5))});
                                }
                            if(cars_closeby.id_r==99 & (cars_closeby.id_r_infront==99 | (cars_closeby.distance_r+cars_closeby.speed_r*4)>(cars_closeby.distance+cars_closeby.speed*4))){
                                goals.push_back({car_s+near_goal, 2*4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_r/5))});
                                goals.push_back({car_s+far_goal, 2*4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_r/5))});
                                }

                            }
                        //when in left lane...
                        if(CURRENT_LANE==0 & cars_closeby.id_r==99 & (cars_closeby.id_r_infront==99 | (cars_closeby.distance_r+cars_closeby.speed_r*4)>(cars_closeby.distance+cars_closeby.speed*4))){
                            goals.push_back({car_s+near_goal, 4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_r/5))});
                            goals.push_back({car_s+far_goal, 4+2, VMAX*1/(1+exp(KEEP_DISTANCE/5-cars_closeby.distance_r/5))});
                            }
                        }
                    }
                //}
            //cout<<"goal lane "<<GOAL_LANE<<" current lane "<<CURRENT_LANE<<" is changing lanes "<<IS_CHANGING_LANES<<endl;
            //cout<<cars_closeby.distance_l+cars_closeby.speed_l*4<<'\t'<<cars_closeby.distance+cars_closeby.speed*4<<'\t'<<cars_closeby.distance_r+cars_closeby.speed_r*4<<endl;

            //check each trajectory for minimum cost
            vector<Trajectory> trajectories;
            double minimum_cost=INFINITY;
            double min_cost_traj=99;
            for(int i=0; i<goals.size();i++){
                    Trajectory each_traj=create_splines_xt_yt(ego_car_data, previous_path_x, previous_path_y, goals[i],  map_waypoints_s,  map_waypoints_x, map_waypoints_y, sensor_fusion );
                trajectories.push_back(each_traj);
                if(each_traj.cost<minimum_cost){
                    minimum_cost=each_traj.cost;
                    min_cost_traj=i;
                }
            }

            for(int i=0; i<goals.size();i++){
                    cout<<goals[i][1]<<" in "<< goals[i][0]-ego_car_data[2]<<"m @ "<<goals[i][2]<<" m/s, Lane-change: " <<trajectories[i].includes_lanechange<<" "<<trajectories[i].cost<<endl;
            }
            //cout<<"Picked "<<min_cost_traj<<endl;


            Trajectory next_traj = trajectories[min_cost_traj];
            if(next_traj.includes_lanechange==1){
                IS_CHANGING_LANES=1;
                //cout<<"about to change_lanes"<<endl;
                GOAL_LANE=round((goals[min_cost_traj][1]-2)/4);
                //cout<<"Goal_lane "<<GOAL_LANE<<endl;
            }
            //GOAL_TIME=next_traj.t_action_complete+TIME;

            if(abs(ego_car_data[3]-(GOAL_LANE*4+2))<0.5){
          	//if(TIME>GOAL_TIME){
                IS_CHANGING_LANES=0;
                //cout<<"overtaking complete"<<endl;
          	}

          	cout<<endl;


          for(int i = 0; i < path_size; i++)
          {
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);
          }


          double how_much_to_add= std::fmin(30.0,(next_traj.t_action_complete/OVERALL_DT)); //maximum of 40  points added to the path


          for(int i = path_size+1; i < how_much_to_add; i++)
          {
              next_x_vals.push_back(next_traj.xt(OVERALL_DT*(i)));
              next_y_vals.push_back(next_traj.yt(OVERALL_DT*(i)));
          }


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
