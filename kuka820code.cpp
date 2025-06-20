
#include <iostream>
#include "windows.h"
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include "ik.h"
#include <vector>
#include <Eigen/Dense> 
//#include"alo.cpp"
//#include"ik.cpp"
//#include"eigen-3.4.0/Eigen/src/Core/DenseBase.h"
//#include "eigen3/Eigen/src/Geometry/"
extern "C" {
#include "extApi.h"
}

using namespace std;
#define PI 3.14




/*void fitness(int x1, int x2, int y1, int y2, int z1, int z2) {
	sqrt((x2 - x1) ^ 2 + (y2 - y1) ^ 2 + (z2 - z1) ^ 2);
}*/
void mulMat(double mat1[][4], double mat2[][4])
{
	double rslt[4][4];

	cout << "Multiplication of given two matrices is:\n";

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			rslt[i][j] = 0;

			for (int k = 0; k < 4; k++) {
				rslt[i][j] += mat1[i][k] * mat2[k][j];
			}

			cout << rslt[i][j] << "\t";
		}

		cout << endl;
	}
}

void inverseKinematics() {
	//defining end-effector parametrs


	double  db = 0.360, d1 = 0.42, d3 = 0.4, dn = 0.2;
	double theta1, theta2, theta3=0, theta4, theta5, theta6,theta7;
	double temp5, temp6, temp7;
	double e1, e2, e3, e4;
	/*int Rd[3][3] = {cos(q1), -sin(q1) * cos(q2),sin(q1) * sin(q2),
				   sin(q1),cos(q1) * cos(q2), -cos(q1) * sin(q2),
				   0,sin(q2),cos(q2) };
	int th[1][3] = {0,0,1};*/
	double xd =0.3;
	double yd =-0.6;// -0.375;//-0.025;
	double zd = 0.2;
	double Od[1][3] = { xd, yd, zd };
	//int Ocd = Od - (d6 * Rd[3][3] * th[1][3]);
	//End - Effector Orientation wrt world frame Rd
	double xc = xd +dn;// *cos(45);// -dn;// *0.7071;
	double yc = yd;// -dn;// *cos(45);// +dn;// +dn;// -dn;// *sin(45));
	double zc = zd-dn;// *cos(45);// +dn;// -(dn * cos(180) * cos(45)); when lower that 0.3 , - in higher than 
	double rc = pow(pow(yc, 2) + pow(xc, 2), 0.5f);
	double s = zc -db;
	double t1 = atan2(xc,yc);
	theta1 = t1 * (180 / PI);
	double D = (pow(rc, 2) + pow(s, 2) - pow(d1, 2) - pow(d3, 2)) / 2 * d1 * d3;
	double t4 = acos(D);
	double t41 = atan2(D, -sqrt(1 - (D * D)));  //acos(D); 
	double theta41 = t41 * (180 / PI);
	double t42 = atan2(D, sqrt(1 - (D * D)));
	double theta42 = t42 * (180 / PI);
	theta4 = t4 * (180 / PI);
	double t2 = atan2(s, rc) - atan2(d3 * sin(theta42), d1 + (d3 * cos(theta42)));
	theta2 = t2 * (180 / PI);
	double co1 = cos(theta1); double co2 = cos(theta2); double co3 = cos(theta4);
	double s1 = sin(theta1); double s2 = sin(theta2); double s3 = sin(theta4);	
	
	
	
	double R11 = co3 * co2;//xyz
	double R12 = s3* co2;
	double R13 = -s2;
	double R21 = -s3 * co1 + co3 * s2 * s1;
	double R22 = co1 * co3 + s1 * s2 * s3;
	double R23 =  co2* s1;
	double R31 = s1 * s3 + co1 * s2 * co3;
	double R32 = -co3 * s1 + co1 * s2 * s3;
	double R33 = co2 * co1;
		/*double R11 = co3 * co2;
	double R12 = -s3 * co1 + co3 * s2 * s1;
	double R13 = s3 * s1 + co3 * s2 * co1;
	double R21 = s3 * co2;
	double R22 = co3 * co1 + s3 * s2 * s1;
	double R23 = -co3 * s1 + s3 * s2 * s1;
	double R31 = -s2;
	double R32 = co2 * s1;
	double R33 = co2 * co1; */
/*double R11 = co2;//xzx
	double R12 = -s2 * co3;
	double R13 = s3 * s2;
	double R21 = s2 * co1;
	double R22 = co1*co2*co3-s3*s1;
	double R23 = -s3*co1*co2-s1*co3;
	double R31 = s1*s2;
	double R32 = s1*co2*co3+s3*co1;
	double R33 = -s1*co2*s3+co1*co3;*/
	
	
	double alph =  acos((R11 + R22 + R33 - 1) / 2);
	double k1 = (R32 - R23) / (2 * sin(alph));
	double k2 = (R13 - R31) / (2 * sin(alph));
	double k3 = (R21 - R12) / (2 * sin(alph));
	e1 = k1* sin(alph / 2);
	e2 = k2* sin(alph / 2);
	e3 =  k3 * sin(alph / 2);
	e4 = cos(alph / 2);
	double cR11 = 1 - 2 * pow(e2, 2) - 2 * pow(e3, 2);
	double cR12 = 2 * (e1 * e2 - e3 * e4);
	double cR13 = 2 * (e1 * e3 + e2 * e4);
	double cR21 = 2 * (e1 * e2 + e3 * e4);
	double cR22 = 1 - 2 * pow(e1, 2) - 2 * pow(e3, 2);
	double cR23 = 2 * (e2 * e3 - e1 * e4);
	double cR31 = 2 * (e1 * e3 - e2 * e4);
	double cR32 = 2 * (e2 * e3 + e1 * e4);
	double cR33= 1 - 2 * pow(e1, 2) - 2 * pow(e2, 2); 
	theta5 = atan2(2 * (e1 * e2 - e4 * e3), e4 * e4 - e3 * e3 + e2 * e2 - e1 * e1);
	theta6 = asin(2 * e2 * e3 + 2 * e4 * e1);
	theta7 = atan2(2 * (e1 * e3 - e4 * e2), e4 * e4 + e3 * e3 - e2 * e2 - e1 * e1);
	
		double res = pow(e1, 2) + pow(e2, 2) + pow(e3, 2) + pow(e4, 2);
		cout << res << endl;
	
	/*theta6 = atan2(-R31, sqrt(R11 * R11 + R21 * R21)) * 180 / PI;
		theta7 = atan2(R32 / cos(theta6), R33 / cos(theta6)) * 180 / PI; 
		theta5 = atan2(R21 / cos(theta6), R11 / cos(theta6)) * 180 / PI; 
*/	 
	/*double t5 = atan2(R23, R13);
	temp5 = t5 * (180 / PI);
	double t6 = acos(R33);  //acos(R33);//atan2(sqrt(R23*R23+R13*R13), R33); 
	temp6 = t6 * (180 / PI);
	temp7 = atan2(R32, -R31) * (180 / PI);
	if (temp5 >= -170 && temp5 <= 170 && temp6 >= -120 && temp6 <= 120 && temp7 >= -175 && temp7 <= 175) {
		theta5 = temp5; 
		theta6 = temp6; 
		theta7 = temp7; 
	} 
	else { 
		 theta5 = atan2(-R23, -R13)*180/PI;
	 	theta6 = -acos(R33)*180/PI;
		theta7 = atan2(-R32, R31)*180/PI;}*/

		/*theta5 = atan2(2 * (e1 * e2 - e4 * e3), e4 * e4 - e3 * e3 + e2 * e2 - e1 * e1);
	theta6 = asin(2 * e2 * e3 + 2 * e4 * e1);
	theta7= atan2(2 * (e1 * e3 - e4 * e2), e4 * e4 +e3 * e3 - e2 * e2 - e1 * e1);
*/
		cout << theta1 <<" "<< theta2 <<" "<< theta3 <<" "<< theta4 <<" "<< theta5 << " "<<theta6 <<" "<< theta7<<endl;
	//cout << "Arm Link: theta 1 ->" << theta1 << " theta 2 ->" << theta2 << " theta 3 ->" << theta3 << "theta 4 ->" << theta4 << endl; 
	//cout << "Wrist Link: theta 5->" << theta5 << " theta 6 ->" << theta6 << " theta 7 ->" << theta7 << endl;
	return;
	
} 
void forwardKinematics(double d1, double d2, double d3, double d4, double d5, double d6, double d7,
	double q1, double q2, double q3, double q4, double q5, double q6, double q7) {
	double c1 = cos(q1); double c2 = cos(q2); double c3 = cos(q3); double c4 = cos(q4); double c5 = cos(q5); double c6 = cos(q6); double c7 = cos(q7);
	double s1 = sin(q1); double s2 = sin(q2); double s3 = sin(q3); double s4 = sin(q4); double s5 = sin(q5); double s6 = sin(q6); double s7 = sin(q7);

	double R11 = c7 * s6 * (s4 * (s1 * s3 - c1 * c2 * c3) - c1 * c4 * s2) * c6 * (c5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) + s5 * (c3 * s1 + c1 * c2 * s3)) + s7 * (s5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) - c5 * (c3 * s1 - c1 * c1 * s3));

	double R12 = c7 * (s5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) - c5 * (c3 * s1 + c1 * c2 * c3)) - s7 * (s6 * (s4 * (s1 * s3 - c1 * c2 * c3) - c1 * c4 * s2) - c6 * (c5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) + s5 * (c3 * s1 + c1 * c2 * c3)));
	double R13 = -c6 * (s4 * (s1 * s3 - c1 * c2 * c3) - c1 * c4 * s2) - s6 * (c5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) + s5 * (c3 * s1 + c1 * c2 * s3));
	double R21 = -(c7 * (s6 * (s4 * (c1 * s3 + c2 * c3 * s1) + c4 * s1 * s2) - c6 * (c5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) + c5 * (c1 * c3 - c2 * s1 * s3)))) - s7 * (s5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) - c5 * (c1 * c3 - c2 * s1 * s3));
	double R22 = -(s7 * (s6 * (s4 * (c1 * s3 + c2 * c3 * s1) + c4 * s1 * s2) - c6 * (c5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) + c5 * (c1 * c3 - c2 * s1 * s3)))) - s7 * (s5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) - c5 * (c1 * c3 - c2 * s1 * s3));
	double R23 = c6 * (s4 * (c1 * s3 + c2 * c3 * s1) + c4 * s1 * s2) + s6 * (c5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) - c5 * (c1 * c3 - c2 * s1 * s3));
	double R33 = c6 * (c2 * c4 - c3 * s2 * s4) - s6 * (c5 * (c2 * s4 + c3 * c4 * s2) - s2 * s3 * s5);
	double Px = d7 * (s4 * (s1 * s3 - c1 * c2 * c3) - c1 * c4 * s2) - d7 * (c6 * (s4 * (s1 * s3 - c1 * c2 * c3) - c1 * c4 * s2) + s6 * (c5 * (c4 * (s1 * s3 - c1 * c2 * c3) + c1 * s2 * s4) + s5 * (c3 * s1 + c1 * c2 * s3)));
	double Py = (d7 * (c6 * (s4 * (c1 * s3 + c2 * c3 * s1) + c4 * s1 * s2) + s6 * (c5 * (c4 * (c1 * s3 + c2 * c3 * s1) - s1 * s2 * s4) + s5 * (c1 * c3 - c2 * s1 * s3)))) * d5 * (s4 * (c1 * s3 + c2 * c3 * s1) + c4 * s1 * s2) + d3 * s1 * s2;
	double Pz = d1 + d5 * (c2 * c4 - c3 * c2 * c4) - (d7 * (s6 * (c5 * (c2 * s4 + c3 * c4 * s2) - s2 * s3 * s5) - s6 * (c2 * c4 - c3 * s2 * s4)) + d3 * c2);
	double T07[4][4] = { R11,R12,R13,Px ,
						 R21,R22,R23,Py ,
						 0,  -1, R33,Pz ,
						  0,   0, 0,  1 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			cout << T07[i][j] << endl;
		}
	}
	return;

}

// simxSetJointTargetPosition(ClientID, lbrJoint1, 0.0, simx_opmode_oneshot_wait);
int main()
{

	//ikCreateEnvironment;
	//int tipHandle, targetHandle; 
	//int joint1Handle, joint2Handle, joint3Handle; 
	//bool ikGetObjectHandle("tip");
	//bool ikGetObjectHandle("target"); 
	//bool ikCreateGroup("env");


	bool VERBOSE = true;
	int ClientID = 0;
	int leftMotorHandle = 0;
	int rightMotorHandle = 0;

	int lbrJoint1 = 0;
	int lbrJoint2 = 0;
	int lbrJoint3 = 0;
	int lbrJoint4 = 0;
	int lbrJoint5 = 0;
	int lbrJoint6 = 0;
	int lbrJoint7 = 0;
	//int lbrJoint8 = 0;

	float minDist = 0;
	float* ptr = &minDist;
	int counter = 0;
	int portNb = 19000;
	bool WORK = true;
	simxFinish(-1);
	ClientID = simxStart((simxChar*)"127.0.0.1", portNb, false, true, 5000, 5);
	Sleep(1000);
	if (ClientID != -1) {
		cout << " connection status to VREP: SUCCESS" << endl;
		simxInt syncho = simxSynchronous(ClientID, 1);
		int start = simxStartSimulation(ClientID, simx_opmode_oneshot_wait);
		int TEST1 = simxGetObjectHandle(ClientID, "Pioneer_p3dx_leftMotor", &leftMotorHandle, simx_opmode_oneshot_wait);
		int TEST2 = simxGetObjectHandle(ClientID, "Pioneer_p3dx_rightMotor", &rightMotorHandle, simx_opmode_oneshot_wait);

		simxSetJointTargetVelocity(ClientID, lbrJoint1, 200, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint2, 500, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint3, 90, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint4, 400, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint5, 300, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint6, 100, simx_opmode_oneshot);
		simxSetJointTargetVelocity(ClientID, lbrJoint7, 200, simx_opmode_oneshot);
		//simxCheckDistance(ClientID, TEST1, TEST2, simxFloat * minimumDistance, simx_opmode_streaming);
		//simxInt simxSetVisionSensorImage(simxInt clientID, simxInt sensorHandle, simxUChar * image, simxInt bufferSize, simxUChar options, simxInt operationMode);
		int target = simxCreateDummy(ClientID, 0.1, nullptr, &lbrJoint7, simx_opmode_blocking);
		//ikCreateEnvironment();
		//int sensor = simxGetObjectHandle(ClientID, "Vision_sensor", nullptr, simx_opmode_oneshot);
	//	simxGetVisionSensorImage(ClientID, sensor, 0, 0, 0, simx_opmode_streaming);
	//	if (simxGetVisionSensorImage(ClientID, sensor, 0, 0, 0, simx_opmode_streaming) == simx_return_ok)
	//		int imageAcquisitionTime = simxGetLastCmdTime(ClientID); 


		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint1", &lbrJoint1, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint2", &lbrJoint2, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint3", &lbrJoint3, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint4", &lbrJoint4, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint5", &lbrJoint5, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint6", &lbrJoint6, simx_opmode_oneshot_wait);
		simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint7", &lbrJoint7, simx_opmode_oneshot_wait);
		//simxGetObjectHandle(ClientID, "LBR_iiwa_14_R820_joint8", &lbrJoint8, simx_opmode_oneshot_wait);

		if (VERBOSE)
		{
			cout << "computer object handle:" << TEST1 << " " << leftMotorHandle << endl;
			cout << "computer object handle" << TEST2 << " " << rightMotorHandle << endl;

		}

		simxSetJointPosition(ClientID, lbrJoint1, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint2, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint3, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint4, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint5, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint6, 0.0, simx_opmode_oneshot);
		simxSetJointPosition(ClientID, lbrJoint7, 0.0, simx_opmode_oneshot);

		cout << "At second Block..." << endl;
		//ikCreateDummy("tip", &lbrJoint1); 
		/*double q = 45 * PI / 180;
		double T1[4][4] = { 1,       0 ,       0 ,    0,
				   0,   cos(q), -sin(q),    0,
				   0,   sin(q),   cos(q),    0,
				   0,        0,        0,    1 };
		double T2[4][4] = { cos(q), -sin(q), 0,    0,
				   sin(q), cos(q),  0,    0,
						0,      0,  1,    0,
						0,      0,  0,    1 };

		mulMat(T1, T2);
*/
		cout << "inverse kinematics" << endl;
		inverseKinematics();

		cout << "end" << endl;
		//simxGetJointMatrix(ClientID, lbrJoint1, , simx_opmode_streaming);

		double q1, q2, q3, q4, q5, q6, q7;
		cout << "Enter 7 angles ";
		//cout << "q0: "; cin >> q1;
		cin >> q1 >> q2 >> q3 >> q4 >> q5 >> q6 >> q7;
		/*/cout << "q1: "; cin >> q1;
		cout << "q2: "; cin >> q2;
		cout << "q3: "; cin >> q3;
		cout << "q4: "; cin >> q4;
		cout << "q5: "; cin >> q5;
		cout << "q6: "; cin >> q6;
		cout << "q7: "; cin >> q7;*/

		cout << "Forward kinematics" << endl;
		double  d1 = 0.360, d3 = 0.42, d5 = 0.4, d7 = 0.2, d2 = 0, d4 = 0, d6 = 0;
		forwardKinematics(d1, d2, d3, d4, d5, d6, d7,
			q1 * (PI / 180), q2 * (PI / 180), q3 * (PI / 180), q4 * (PI / 180), q5 * (PI / 180), q6 * (PI / 180), q7 * (PI / 180));


		//simxSetJointTargetPosition(ClientID, lbrJoint1, q0 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint1, q1 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint2, q2 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint3, q3 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint4, q4 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint5, q5 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint6, q6 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint7, q7 * (PI / 180), simx_opmode_oneshot_wait);

		//simxGetJointPosition(ClientID, lbrJoint7, 0, simx_opmode_streaming);
		///////////peredelat have to make a choice as Game.cpp
		cout << " If you want to change a join to cajmge choose:";


		while (true) {

			if (GetAsyncKeyState((unsigned short)'A') & 0x8000) {
				q1 += 15;
			}
			if (GetAsyncKeyState((unsigned short)'D') & 0x8000) {
				q1 -= 15;
			}
			if (GetAsyncKeyState((unsigned short)'W') & 0x8000) {
				q2 += 15;
			}
			if (GetAsyncKeyState((unsigned short)'S') & 0x8000) {
				q2 -= 15;
			
			}if (GetAsyncKeyState((unsigned short)'J') & 0x8000) {
				q4 += 15;
			}
			if (GetAsyncKeyState((unsigned short)'L') & 0x8000) {
				q4 -= 15;
			}
			if (GetAsyncKeyState((unsigned short)'R') & 0x8000) {
				q1 = 0;
				q2 = 0;
				q3 = 0;
				q4 = 0;
				q5 = 0;
				q6 = 0;
				q7 = 0;
			}
		simxSetJointTargetPosition(ClientID, lbrJoint1, q1 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint2, q2 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint3, q3 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint4, q4 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint5, q5 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint6, q6 * (PI / 180), simx_opmode_oneshot_wait);
		simxSetJointTargetPosition(ClientID, lbrJoint7, q7 * (PI / 180), simx_opmode_oneshot_wait);
	}

		////////////
		cout << simxCheckDistance(ClientID, leftMotorHandle, rightMotorHandle, &minDist, simx_opmode_streaming) << endl;

		simxGetJointPosition(ClientID, leftMotorHandle, 0, simx_opmode_streaming);
		simxGetJointPosition(ClientID, rightMotorHandle, 0, simx_opmode_streaming);
		
		//bool ikCreateEnviroment(nullptr);
		//simxGetJointMatrix(ClientID,leftMotorHandle, simxFloat * matrix, simx_opmode_streaming);
		 //simGetObject(const simChar * objectPath, simInt index, simInt proxy, simInt options);
		//simxGetObjectGroupData( ClientID, sim_appobj_object_type, sim_object_joint_type, nullptr, nullptr, simxInt * intDataCount, simxInt * *intData, simxInt * floatDataCount, simxFloat * *floatData, simxInt * stringDataCount, simxChar * *stringData, simx_opmode_blocking)
	}
	else {
		cout << "Connection status to VREP: FAILED" << endl;
	}
	simxFinish(ClientID);
	return ClientID;
}
