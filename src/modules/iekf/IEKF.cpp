#include "IEKF.hpp"

static const float mag_inclination = 1.0f;
static const float mag_declination = 0;

IEKF::IEKF() :
	_nh(), // node handle
	_sub_gyro(_nh.subscribe("sensor_gyro", 0, &IEKF::callback_gyro, this)),
	_sub_accel(_nh.subscribe("sensor_accel", 0, &IEKF::callback_accel, this)),
	_sub_mag(_nh.subscribe("sensor_mag", 0, &IEKF::callback_mag, this)),
	_sub_baro(_nh.subscribe("sensor_baro", 0, &IEKF::callback_baro, this)),
	_sub_gps(_nh.subscribe("vehicle_gps_position", 0, &IEKF::callback_gps, this)),
	_pub_attitude(_nh.advertise<vehicle_attitude_s>("vehicle_attitude", 0)),
	_pub_local_position(_nh.advertise<vehicle_local_position_s>("vehicle_local_position", 0)),
	_pub_global_position(_nh.advertise<vehicle_global_position_s>("vehicle_global_position", 0)),
	_pub_control_state(_nh.advertise<control_state_s>("control_state", 0)),
	_x(),
	_P(),
	_u(),
	_g_n(0, 0, -9.8),
	_B_n(),
	_map_ref()
{
	// start with 0 quaternion
	_x(X::q_nb_0) = 1;
	_x(X::q_nb_1) = 0;
	_x(X::q_nb_2) = 0;
	_x(X::q_nb_3) = 0;

	// start with 1 accel scale
	_x(X::accel_scale) = 1;

	// initialize covariance
	_P(Xe::rot_n, Xe::rot_n) = 1e-2f;
	_P(Xe::rot_e, Xe::rot_e) = 1e-2f;
	_P(Xe::rot_d, Xe::rot_d) = 1e-2f;
	_P(Xe::vel_n, Xe::vel_n) = 1;
	_P(Xe::vel_e, Xe::vel_e) = 1;
	_P(Xe::vel_d, Xe::vel_d) = 1;
	_P(Xe::gyro_bias_n, Xe::gyro_bias_n) = 1e-2f;
	_P(Xe::gyro_bias_e, Xe::gyro_bias_e) = 1e-2f;
	_P(Xe::gyro_bias_d, Xe::gyro_bias_d) = 1e-2f;
	_P(Xe::accel_scale, Xe::accel_scale) = 1;
	_P(Xe::pos_n, Xe::pos_n) = 1;
	_P(Xe::pos_e, Xe::pos_e) = 1;
	_P(Xe::pos_d, Xe::pos_d) = 1;
	_P(Xe::terrain_alt, Xe::terrain_alt) = 1;
	_P(Xe::baro_bias, Xe::baro_bias) = 1;

	// initial magnetic field guess
	_B_n = Vector3f(0.21523, 0.00771, -0.42741);
}

Vector<float, X::n> IEKF::dynamics(const Vector<float, X::n> &x, const Vector<float, U::n> &u)
{
	Quaternion<float> q_nb(x(X::q_nb_0), x(X::q_nb_1), x(X::q_nb_2), x(X::q_nb_3));
	Vector3<float> a_b(_u(U::accel_bx), _u(U::accel_by), _u(U::accel_bz));
	Vector3<float> as_n = q_nb.conjugate(a_b / _x(X::accel_scale)) - _g_n;
	Vector3<float> gyro_bias_b(_x(X::gyro_bias_bx), _x(X::gyro_bias_by), _x(X::gyro_bias_bz));
	Vector3<float> omega_nb_b(_u(U::omega_nb_bx), _u(U::omega_nb_by), _u(U::omega_nb_bz));
	Quaternion<float> dq_nb = q_nb.derivative(omega_nb_b - gyro_bias_b);

	Vector<float, X::n> dx;
	dx.setZero();
	dx(X::q_nb_0) = dq_nb(0);
	dx(X::q_nb_1) = dq_nb(1);
	dx(X::q_nb_2) = dq_nb(2);
	dx(X::q_nb_3) = dq_nb(3);
	dx(X::vel_n) = as_n(0);
	dx(X::vel_e) = as_n(1);
	dx(X::vel_d) = as_n(2);
	dx(X::gyro_bias_bx) = 0;
	dx(X::gyro_bias_by) = 0;
	dx(X::gyro_bias_bz) = 0;
	dx(X::accel_scale) = 0;
	dx(X::pos_n) = x(X::vel_n);
	dx(X::pos_e) = x(X::vel_e);
	dx(X::pos_d) = x(X::vel_d);
	dx(X::terrain_alt) = 0;
	dx(X::baro_bias) = 0;
	return dx;
}

void IEKF::callback_gyro(const sensor_gyro_s *msg)
{
	//ROS_INFO("gyro callback %10.4f %10.4f %10.4f",
	//double(msg->x), double(msg->y), double(msg->x));
	_u(U::omega_nb_bx) = msg->x;
	_u(U::omega_nb_by) = msg->y;
	_u(U::omega_nb_bz) = msg->z;

	// predict driven by gyro callback
	if (msg->integral_dt > 0) {
		predict(msg->integral_dt / 1.0e6f);
	};
}

void IEKF::callback_accel(const sensor_accel_s *msg)
{
	//ROS_INFO("accel callback %10.4f %10.4f %10.4f",
	//double(msg->x), double(msg->y), double(msg->x));
	_u(U::accel_bx) = msg->x;
	_u(U::accel_by) = msg->y;
	_u(U::accel_bz) = msg->z;

	// calculate residual
	Quaternion<float> q_nb(_x(X::q_nb_0), _x(X::q_nb_1),
			       _x(X::q_nb_2), _x(X::q_nb_3));
	Vector3f y_b(msg->x, msg->y, msg->z);
	Vector3f r = q_nb.conjugate(y_b / _x(X::accel_scale)) - _g_n;

	// define R
	Matrix<float, Y_accel::n, Y_accel::n> R;
	R(Y_accel::accel_bx, Y_accel::accel_bx) = 100.0;
	R(Y_accel::accel_by, Y_accel::accel_by) = 100.0;
	R(Y_accel::accel_bz, Y_accel::accel_bz) = 100.0;

	// define H
	Matrix<float, Y_accel::n, Xe::n> H;
	Matrix3f tmp = _g_n.hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			H(Y_accel::accel_bx + i, Xe::rot_n + j) = tmp(i, j);
		}
	}

	bool fault = correct<Y_accel::n>(r, H, R);

	if (fault) {
		ROS_WARN("accel fault");
	}
}

void IEKF::callback_mag(const sensor_mag_s *msg)
{
	//ROS_INFO("mag callback %10.4f %10.4f %10.4f",
	//double(msg->x), double(msg->y), double(msg->x));

	// calculate residual
	Quaternion<float> q_nb(_x(X::q_nb_0), _x(X::q_nb_1),
			       _x(X::q_nb_2), _x(X::q_nb_3));
	Vector3<float> y_b = Vector3<float>(msg->x, msg->y, msg->z).unit();
	Vector3<float> B_n = _B_n.unit();
	Vector3<float> r = q_nb.conjugate(y_b) - B_n;

	// define R
	Matrix<float, Y_mag::n, Y_mag::n> R;
	R(Y_mag::mag_n, Y_mag::mag_n) = 100;
	R(Y_mag::mag_e, Y_mag::mag_e) = 100;
	R(Y_mag::mag_d, Y_mag::mag_d) = 10000.0; // prevents mag from correcting roll/pitch

	// define H
	Matrix<float, Y_mag::n, Xe::n> H;
	Matrix3f tmp = B_n.hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			H(Y_mag::mag_n + i, Xe::rot_n + j) = tmp(i, j);
		}
	}

	bool fault = correct<Y_mag::n>(r, H, R);

	if (fault) {
		ROS_WARN("mag fault");
	}
}

void IEKF::callback_baro(const sensor_baro_s *msg)
{
	//ROS_INFO("baro callback %10.4f", double(msg->altitude));

	// calculate residual
	Vector<float, Y_baro::n> y;
	y(Y_baro::asl) = msg->altitude;
	Vector<float, Y_baro::n> yh;
	yh(Y_baro::asl)	= -_x(X::pos_d) + _x(X::baro_bias);
	Vector<float, Y_baro::n> r = y - yh;

	// define R
	Matrix<float, Y_baro::n, Y_baro::n> R;
	R(Y_baro::asl, Y_baro::asl) = 1.0;

	// define H
	Matrix<float, Y_baro::n, Xe::n> H;
	H(Y_baro::asl, Xe::pos_d) = -1;
	H(Y_baro::asl, Xe::baro_bias) = 1;

	bool fault = correct<Y_baro::n>(r, H, R);

	if (fault) {
		ROS_WARN("baro fault");
	}
}

void IEKF::callback_gps(const vehicle_gps_position_s *msg)
{
	// check for good gps signal
	if (msg->satellites_used < 6 || msg->fix_type < 3) {
		return;
	}

	// init global reference
	if (!_map_ref.init_done) {
		ROS_INFO("gps map ref init", double(msg->lat * 1e-7), double(msg->lon * 1e-7));
		map_projection_init(&_map_ref,
				    msg->lat * 1e-7, msg->lon * 1e-7);
	}

	// calculate residual
	float pos_n = 0;
	float pos_e = 0;
	map_projection_project(&_map_ref, msg->lat * 1e-7, msg->lon * 1e-7, &pos_n, &pos_e);

	Vector<float, Y_gps::n> y;
	y(Y_gps::pos_n) = pos_n;
	y(Y_gps::pos_e) = pos_e;
	y(Y_gps::pos_d) = msg->alt * 1e-3;
	y(Y_gps::vel_n) = msg->vel_n_m_s;
	y(Y_gps::vel_e) = msg->vel_e_m_s;
	y(Y_gps::vel_d) = msg->vel_d_m_s;

	Vector<float, Y_gps::n> yh;
	yh(Y_gps::pos_n) = _x(X::pos_n);
	yh(Y_gps::pos_e) = _x(X::pos_e);
	yh(Y_gps::pos_d) = _x(X::pos_d);
	yh(Y_gps::vel_n) = _x(X::vel_n);
	yh(Y_gps::vel_e) = _x(X::vel_e);
	yh(Y_gps::vel_d) = _x(X::vel_d);

	Vector<float, Y_gps::n> r = y - yh;

	// define R
	Matrix<float, Y_gps::n, Y_gps::n> R;
	R(Y_gps::pos_n, Y_gps::pos_n) = 1;
	R(Y_gps::pos_e, Y_gps::pos_e) = 1;
	R(Y_gps::pos_d, Y_gps::pos_d) = 1000000;
	R(Y_gps::vel_n, Y_gps::vel_n) = 1;
	R(Y_gps::vel_e, Y_gps::vel_e) = 1;
	R(Y_gps::vel_d, Y_gps::vel_d) = 1;

	// define H
	Matrix<float, Y_gps::n, Xe::n> H;
	H(Y_gps::pos_n, Xe::pos_n) = 1;
	H(Y_gps::pos_e, Xe::pos_e) = 1;
	H(Y_gps::pos_d, Xe::pos_d) = 1;
	H(Y_gps::vel_n, Xe::vel_n) = 1;
	H(Y_gps::vel_e, Xe::vel_e) = 1;
	H(Y_gps::vel_d, Xe::vel_d) = 1;

	bool fault = correct<Y_gps::n>(r, H, R);

	if (fault) {
		ROS_WARN("gps fault");
		r.print();
	}
}

void IEKF::predict(float dt)
{
	// define process noise matrix
	Matrix<float, Xe::n, Xe::n> Q;
	Q(Xe::rot_n, Xe::rot_n) = 1e-1f;
	Q(Xe::rot_e, Xe::rot_e) = 1e-1f;
	Q(Xe::rot_d, Xe::rot_d) = 1e-1f;
	Q(Xe::vel_n, Xe::vel_n) = 1e-2f;
	Q(Xe::vel_e, Xe::vel_e) = 1e-2f;
	Q(Xe::vel_d, Xe::vel_d) = 1e-2f;
	Q(Xe::gyro_bias_n, Xe::gyro_bias_n) = 5e-3f;
	Q(Xe::gyro_bias_e, Xe::gyro_bias_e) = 5e-3f;
	Q(Xe::gyro_bias_d, Xe::gyro_bias_d) = 5e-3f;
	Q(Xe::accel_scale, Xe::accel_scale) = 1e-2f;
	Q(Xe::pos_n, Xe::pos_n) = 1e-2f;
	Q(Xe::pos_e, Xe::pos_e) = 1e-2f;
	Q(Xe::pos_d, Xe::pos_d) = 1e-2f;
	Q(Xe::terrain_alt, Xe::terrain_alt) = 1e-1f;
	Q(Xe::baro_bias, Xe::baro_bias) = 1e-1f;

	// define A matrix
	Matrix<float, Xe::n, Xe::n> A;
	A.setZero();

	// derivative of rotation error is -0.5 * gyro bias
	A(Xe::rot_n, Xe::Xe::gyro_bias_n) = -0.5;
	A(Xe::rot_e, Xe::Xe::gyro_bias_e) = -0.5;
	A(Xe::rot_d, Xe::Xe::gyro_bias_d) = -0.5;

	// derivative of velocity
	Quaternion<float> q_nb(
		_x(X::q_nb_0), _x(X::q_nb_1),
		_x(X::q_nb_2), _x(X::q_nb_3));
	Vector3<float> a_b(_u(U::accel_bx), _u(U::accel_by), _u(U::accel_bz));
	Vector3<float> J_a_n = q_nb.conjugate(a_b / _x(X::accel_scale));
	Matrix<float, 3, 3> a_tmp = -J_a_n.hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			A(Xe::vel_n + i, Xe::rot_n + j) = a_tmp(i, j);
		}

		A(Xe::vel_n + i, Xe::accel_scale) = -J_a_n(i);
	}

	// derivative of gyro bias
	Vector3<float> omega_nb_b(
		_u(U::omega_nb_bx), _u(U::omega_nb_by), _u(U::omega_nb_bz));
	Vector3<float> gyro_bias_b(
		_x(X::gyro_bias_bx), _x(X::gyro_bias_by), _x(X::gyro_bias_bz));
	Vector3<float> J_omega_n = q_nb.conjugate(omega_nb_b - gyro_bias_b);
	Matrix<float, 3, 3> g_tmp = J_omega_n.hat();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			A(Xe::vel_n + i, Xe::rot_n + j) = g_tmp(i, j);
		}

		A(Xe::vel_n + i, Xe::accel_scale) = -J_a_n(i);
	}

	// derivative of position is velocity
	A(Xe::pos_n, Xe::vel_n) = 1;
	A(Xe::pos_e, Xe::vel_e) = 1;
	A(Xe::pos_d, Xe::vel_d) = 1;

	// propgate state using euler integration
	Vector<float, X::n> dx = dynamics(_x, _u) * dt;
	_x += dx;

	// propgate covariance using euler integration
	Matrix<float, Xe::n, Xe::n> dP = (A * _P + _P * A.T() + Q) * dt;
	_P += dP;

	ros::Time now = ros::Time::now();

	// publish attitude
	{
		vehicle_attitude_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.q[0] = _x(X::q_nb_0);
		msg.q[1] = _x(X::q_nb_1);
		msg.q[2] = _x(X::q_nb_2);
		msg.q[3] = _x(X::q_nb_3);
		msg.rollspeed = _u(U::omega_nb_bx);
		msg.pitchspeed = _u(U::omega_nb_by);
		msg.yawspeed = _u(U::omega_nb_bz);
		_pub_attitude.publish(msg);
	}

	// publish local position
	{
		vehicle_local_position_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.xy_valid = true;
		msg.z_valid = true;
		msg.v_xy_valid = true;
		msg.v_z_valid = true;
		msg.x = _x(X::pos_n);
		msg.y = _x(X::pos_e);
		msg.z = _x(X::pos_d);
		msg.vx = _x(X::vel_n);
		msg.vy = _x(X::vel_e);
		msg.vz = _x(X::vel_d);
		_pub_local_position.publish(msg);
	}

	// publish global position
	{
		double lat = 0;
		double lon = 0;
		double alt = -_x(X::pos_d);
		Euler<float> euler_nb = Quaternion<float>(_x(X::q_nb_0),
					_x(X::q_nb_1), _x(X::q_nb_2), _x(X::q_nb_3));
		map_projection_reproject(&_map_ref, _x(X::pos_n), _x(X::pos_e), &lat, &lon);
		vehicle_global_position_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.time_utc_usec = 0; // TODO
		msg.lat = lat;
		msg.lon = lon;
		msg.alt = alt;
		msg.vel_n = _x(X::vel_n);
		msg.vel_e = _x(X::vel_e);
		msg.vel_d = _x(X::vel_d);
		msg.yaw = euler_nb(2);
		msg.eph = sqrt(_P(Xe::pos_n, Xe::pos_n) + _P(Xe::pos_e, Xe::pos_e));
		msg.epv = _P(Xe::pos_d, Xe::pos_d);
		msg.terrain_alt = _x(X::terrain_alt);
		msg.terrain_alt_valid = true;
		msg.dead_reckoning = false;
		msg.pressure_alt = alt; // TODO
		_pub_global_position.publish(msg);
	}

	// publish control state
	{
		control_state_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.x_acc = 0;
		msg.y_acc = 0;
		msg.z_acc = 0;
		msg.x_vel = _x(X::vel_n);
		msg.y_vel = _x(X::vel_e);
		msg.z_vel = _x(X::vel_d);
		msg.x_pos = _x(X::pos_n);
		msg.y_pos = _x(X::pos_e);
		msg.z_pos = _x(X::pos_d);
		msg.airspeed = 0;
		msg.airspeed_valid = false;
		msg.vel_variance[0] = _P(Xe::vel_n, Xe::vel_n);
		msg.vel_variance[1] = _P(Xe::vel_e, Xe::vel_e);
		msg.vel_variance[2] = _P(Xe::vel_d, Xe::vel_d);
		msg.pos_variance[0] = _P(Xe::pos_n, Xe::pos_n);
		msg.pos_variance[1] = _P(Xe::pos_e, Xe::pos_e);
		msg.pos_variance[2] = _P(Xe::pos_d, Xe::pos_d);
		msg.q[0] = _x(X::q_nb_0);
		msg.q[1] = _x(X::q_nb_1);
		msg.q[2] = _x(X::q_nb_2);
		msg.q[3] = _x(X::q_nb_3);
		msg.delta_q_reset[0] = 0;
		msg.delta_q_reset[1] = 0;
		msg.delta_q_reset[2] = 0;
		msg.delta_q_reset[3] = 0;
		msg.quat_reset_counter = 0;
		msg.roll_rate = _u(U::omega_nb_bx);
		msg.pitch_rate = _u(U::omega_nb_by);
		msg.yaw_rate = _u(U::omega_nb_bz);
		msg.horz_acc_mag = 0;
		_pub_control_state.publish(msg);
	}
}

void IEKF::applyErrorCorrection(Vector<float, Xe::n> d_xe)
{
	Quaternion<float> q_nb(_x(X::q_nb_0), _x(X::q_nb_1), _x(X::q_nb_2), _x(X::q_nb_3));
	Quaternion<float> d_q_nb = Quaternion<float>(0,
				   d_xe(Xe::rot_n), d_xe(Xe::rot_e), d_xe(Xe::rot_d)) * q_nb;
	//ROS_INFO("d_q_nb");
	//d_q_nb.print();
	Vector3<float> d_gyro_bias_b = q_nb.conjugate_inversed(
					       Vector3<float>(d_xe(Xe::gyro_bias_n),
							       d_xe(Xe::gyro_bias_e),
							       d_xe(Xe::gyro_bias_d)));

	// linear term correction is the same
	// as the error correction
	_x(X::q_nb_0) += d_q_nb(0);
	_x(X::q_nb_1) += d_q_nb(1);
	_x(X::q_nb_2) += d_q_nb(2);
	_x(X::q_nb_3) += d_q_nb(3);
	_x(X::gyro_bias_bx) += d_gyro_bias_b(0);
	_x(X::gyro_bias_by) += d_gyro_bias_b(1);
	_x(X::gyro_bias_bz) += d_gyro_bias_b(2);
	_x(X::vel_n) += d_xe(Xe::vel_n);
	_x(X::vel_e) += d_xe(Xe::vel_e);
	_x(X::vel_d) += d_xe(Xe::vel_d);
	_x(X::pos_n) += d_xe(Xe::pos_n);
	_x(X::pos_e) += d_xe(Xe::pos_e);
	_x(X::pos_d) += d_xe(Xe::pos_d);
	_x(X::terrain_alt) += d_xe(Xe::terrain_alt);
	_x(X::baro_bias) += d_xe(Xe::baro_bias);
}
