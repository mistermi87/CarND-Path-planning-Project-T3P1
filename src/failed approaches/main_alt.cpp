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
#include "spline.h"
#include <typeinfo>

using namespace std;

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

//-----------------Global variables
int TIME =0;
vector<double> SIZE_OF_CARS_F {2, 1.2}; //length and width of car in frenet coordinates
double OVERALL_DT =0.02; //timestep between pathpoints
double A_MAX=3.5; //max. Acceleration m/s^2
double PREVIOUS_RECOG=10; //how many of the previous path's timepoints are taken into consideration for the new path
double MPH_IN_MS = 0.44704; //M/h in m/s
double VMAX=50*0.44704; //Vmax in m/s
double V_GOAL=45*0.44704; //V_goal in m/s
double S_FULL=6945.554;// length of the track

int LAST_WAYPOINT;
tk::spline XS;
tk::spline YS;


//Weights
double WEIGHT_DIST=1;   //How far does each trajectory take the car in 4s
double WEIGHT_SPEED=1;  //Does the car go over Vmax? is it going exceptionally slow?
double WEIGHT_ACCL=1;   //Dos the trajectory create high accelerations?
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

vector<double> getXY_fine(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y, const vector<double> &maps_dx, const vector<double> &maps_dy)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int prev2_wp=((prev_wp-1)+maps_s.size())%maps_s.size();

	int next_wp = (prev_wp+1)%maps_x.size();
	int next2_wp= (next_wp+1)%maps_x.size();

	double seg_s = std::fmod(((s-maps_s[prev_wp])+S_FULL),S_FULL);

	double s_prev2_wp= -std::fmod(((maps_s[prev_wp]-maps_s[prev2_wp])+S_FULL),S_FULL);
	double s_prev_wp=0;
    double s_next_wp = std::fmod(((maps_s[next_wp]-maps_s[prev_wp])+S_FULL),S_FULL);
	double s_next2_wp = std::fmod(((maps_s[next2_wp]-maps_s[prev_wp])+S_FULL),S_FULL);

	cout<<prev2_wp<<"\t "<<prev_wp<<"\t "<<next_wp<<" \t"<<next2_wp<<endl;
	cout<<s_prev2_wp<<"\t "<<s_prev_wp<<"\t "<<s_next_wp<<"\t "<<s_next2_wp<<endl;

	if (prev_wp-LAST_WAYPOINT!=0){
        XS.set_points({s_prev2_wp,s_prev_wp,s_next_wp,s_next2_wp},{maps_x[prev2_wp],maps_x[prev_wp],maps_x[next_wp],maps_x[next_wp]});
        YS.set_points({s_prev2_wp,s_prev_wp,s_next_wp,s_next2_wp},{maps_y[prev2_wp],maps_y[prev_wp],maps_y[next_wp],maps_y[next_wp]});
        LAST_WAYPOINT=prev_wp;
	}



    double x_s=XS(seg_s);
    double d_x_int=maps_dx[prev_wp]+seg_s/(s_next_wp-s_prev_wp)*(maps_dx[next_wp]-maps_dx[prev_wp]);
    double x=x_s+d*maps_dx[prev_wp];//d_x_int;


    double y_s=YS(seg_s);
    double d_y_int=maps_dy[prev_wp]+seg_s/(s_next_wp-s_prev_wp)*(maps_dy[next_wp]-maps_dy[prev_wp]);
    double y=y_s+d*maps_dy[prev_wp];//d_y_int;

	return {x,y};

}

vector<double> getXY_superfine(double s, double dc, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y, const vector<double> &maps_dx, const vector<double> &maps_dy)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();
	//cout<<"Waypoints"<<prev_wp<<" "<<wp2<<endl;

    //----using cubic spline to calculate x(s) and y(s) of road
    double seg_s = std::fmod(((s-maps_s[prev_wp])+S_FULL),S_FULL);
    double s_dist = std::fmod(((maps_s[wp2]-maps_s[prev_wp])+S_FULL),S_FULL);
    double s_rel=seg_s/s_dist;

    //cout<<"S: "<<seg_s<<" "<<s_wp2<<endl;

    //x(s)
    double x_0=maps_x[prev_wp];
    double x_1=maps_x[wp2];

    double dy_0=maps_dy[prev_wp]; //dy : normal of x(s) --> gradient: -dy*n=Dx (scaler of the gradient)
    double dy_1=maps_dy[wp2];

    double a=x_0;
    double b=-dy_0;
    double c=2*dy_0+dy_1+3*(x_1-x_0);
    double d = 2*(x_0-x_1)-dy_1-dy_0;


    double x_s=a+b*s_rel+c*s_rel*s_rel+d*s_rel*s_rel*s_rel;

    //double d_y_int=-(b+2*c*s_rel+3*d*s_rel*s_rel);

    double d_x_int=maps_dx[prev_wp]+s_rel*(maps_dx[wp2]-maps_dx[prev_wp]);



    cout<<"X: "<<x_0<<" "<<x_1<<endl;
    cout<<"dy_01: "<<dy_0<<" "<<dy_1<<endl;
    cout<<" a: "<<a<<" b: "<<b<<" c: "<<c<<" d: "<<d<<endl;
    //cout<<"x_s: "<<x_s<<" d_x_int: "<<d_x_int<<" x: "<<x<<endl;

    //y(s)
    double y_0=maps_y[prev_wp];
    double y_1=maps_y[wp2];

    double dx_0=maps_dx[prev_wp]; //dy : normal of x(s) --> gradient: -dy*n=Dx (scaler of the gradient)
    double dx_1=maps_dx[wp2];

    //n= (2.0/s_wp2)*(y_1-y_0)/(dx_0+dx_1);

    a=y_0;
    b=dx_0;
    c=3*(y_1-y_0)-2*dx_0-dx_1;
    d=dx_1+dx_0-2*y_1+2*y_0;

    double y_s=a+b*s_rel+c*s_rel*s_rel+d*s_rel*s_rel*s_rel;
    double d_y_int=maps_dy[prev_wp]+s_rel*(maps_dy[wp2]-maps_dy[prev_wp]);

    //double d_x_int=(b+2*c*s_rel+3*d*s_rel*s_rel);

    double x=x_s+dc*d_x_int;
    double y=y_s+dc*d_y_int;

    cout<<"Y: "<<y_0<<" "<<y_1<<endl;
    cout<<"dx_01: "<<dx_0<<" "<<dx_1<<endl;
    cout<<" a: "<<a<<" b: "<<b<<" c: "<<c<<" d: "<<d<<endl;
    //cout<<"y_s: "<<y_s<<" d_y_int: "<<d_y_int<<" y: "<<y<<endl;

    cout<<"s_rel: " <<s_rel<<endl;

	return {x,y};

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

          	vector<double> ego_car_data{car_x,car_y,car_s,car_d,car_yaw,car_speed*0.44704}; //<-- speed conversion to m/s;
            cout<<"--- --- ---"<<endl;
            cout<<" X: "<<ego_car_data[0]<<" Y: "<<ego_car_data[1]<<" S: "<<ego_car_data[2]<<" D: "<<ego_car_data[3]<<" Yaw: "<<ego_car_data[4]<<" v: "<<ego_car_data[5]<<endl;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            int path_size=previous_path_x.size();

          	if(path_size<(PREVIOUS_RECOG)){
                for(int i=0; i<PREVIOUS_RECOG-path_size;i++){
                    previous_path_x.insert(previous_path_x.begin(),car_x);
                    previous_path_y.insert(previous_path_y.begin(),car_y);
                }
            }

            for(int i = 0; i < 0; i++)
            {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }
            //cout<<endl;

            double dist_inc = 0.5;

            double s=car_s;
            double d=6;
            vector<double> XY_from_quadint;

            for(int i = 0; i < (50-0); i++)
            {
                XY_from_quadint=getXY_superfine(s+(i+0)*dist_inc, d, map_waypoints_s, map_waypoints_x, map_waypoints_y, map_waypoints_dx, map_waypoints_dy);
                next_x_vals.push_back(XY_from_quadint[0]);
                next_y_vals.push_back(XY_from_quadint[1]);
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
