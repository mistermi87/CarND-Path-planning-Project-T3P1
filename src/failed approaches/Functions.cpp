#include <math.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "JMT.cpp"
#include "spline.h"


vector<tk::spline> create_splines_st_dt(vector<double> car_data, vector<double> goal, vector<double> map_waypoints_x, vector<double> map_waypoints_y ){

    double goal_s=goal[0];
    double goal d=goal[1];
    double goal_velocity=goal[2];

    vector<double> loc_last_sec=getFrenet(car_data[0]-cos(deg2rad(car_data[4]))*car_data[5], car_data[1]+sin(deg2rad(car_data[4]))*car_data[5], car_data[4], map_waypoints_x, map_waypoints_y)

    double d_m_1=loc_last_sec[1]
    double d_0=car_data[3]
    double d_1=goal_d;
    double d_2=goal_d;


    double s_m_1=loc_last_sec[0]
    double s_0=car_data[2]
    double s_1=goal_s
    double s_2=goal_s+0.1*goal_velocity

    double t_m_1=-1;
    double t_0=0;
    double t_1=goal_s/goal_velocity
    double t_2=t_1+0.1;

    tk::spline st;
    st.set_points({t_m_1,t_0,t_1,t_2},{s_m_1,s_0,s_1,s_2});

    tk::spline dt;
    dt.set_points({t_m_1,t_0,t_1,t_2},{d_m_1,d_0,d_1,d_2});

    return {st,dt}
}


vector<tk::spline> create_splines_xt_yt(vector<double> car_data, vector<double> goal, vector<double> map_waypoints_s, vector<double> map_waypoints_x, vector<double> map_waypoints_y ){



    vector<double> goalXY= getXY(goal[0], goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y)

    double goal_velocity=goal[2];
    double goal_dt = 0.1;

    vector<double> goal_2_XY= getXY(goal[0]+goal_dt*goal_velocity, goal[1], map_waypoints_s, map_waypoints_x, map_waypoints_y)

    double x_m_1=car_data[0]-cos(deg2rad(car_data[4]))*car_data[5];
    double x_0=car_data[0]
    double x_1=goalXY[0];
    double d_2=goal_2_XY[0];


    double y_m_1=car_data[1]-sin(deg2rad(car_data[4]))*car_data[5];
    double y_0=car_data[1]
    double y_1=goalXY[1]
    double y_2=goal_2_XY[1]

    double t_m_1=-1;
    double t_0=0;
    double t_1=distance(car_data[0],car_data[1],goalXY[0],goalXY[1])/goal_velocity;
    double t_2=t_1+goal_dt;

    tk::spline xt;
    st.set_points({t_m_1,t_0,t_1,t_2},{x_m_1,x_0,x_1,x_2});

    tk::spline yt;
    dt.set_points({t_m_1,t_0,t_1,t_2},{y_m_1,y_0,y_1,y_2});

    return {xt,yt}
}
