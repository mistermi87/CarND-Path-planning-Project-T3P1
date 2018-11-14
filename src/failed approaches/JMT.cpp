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

    cout<<endl;
	for(int i = 0; i < 6; i++)
	{
	    cout<<result[i]<<" ";
	}
	cout<<endl;

    return result;
}
