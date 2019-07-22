#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "kengine.hpp"

using namespace std;
SDL_Event theEvent;

int mausx;
int mausy;
int runs;

int t1;
int t2;
int t3;
int t4;

double mx;
double my;
double lastmx;
double lastmy;
double a;
double b;
double c;
double d;
double x;
double p1x = 200;
double p1y = 250;
double p2x = 523;
double p2y = 433;
double p1p2vx = (p2y - p1y)/5;
double p1p2vy = -((p2x - p1x)/5);
double s;
double s1;
double s2;

KE_Material *regBall;
KE_Material *noGravBall;
KE_Ball *tb[12];
KE_Edge *testEdge;
KE_Angle *ta[2];

int main(int argc, char* args[]){

	KA_InitAll("blabla", 800, 600, 800, 600, 0);

	KG_LoadTexturesWithTxdFile("images/textures.txt");
	KG_LoadTexturesWithTxdFile("images/textures2.txt");

    regBall = KP_CreateMaterial(0.8, 1);
    noGravBall = KP_CreateMaterial(1, 1, 0xffffffff, 0, 0);

    tb[9] = KP_CreateBall(400, 315, M_PI_2, regBall, 14, 10);

	for (int i = 0; i < 9; i++){
		tb[i] = KP_CreateBall(200+(i%3)*30, 270+(i/3)*30+i, 0, regBall, 13, 10);
	}

    tb[10] = KP_CreateBall(100, 400, 0, noGravBall, 13, 10000000);
    tb[11] = KP_CreateBall(400, 500, 0, noGravBall, 13, 10000000);

	for (int i = 0; i < 12; i++){
		KG_AddBallBindedTexture(7, 1, tb[i], 13);
	}

    testEdge = KP_CreateEdge(tb[11], tb[10], noGravBall, 1);
    ta[0] = KP_CreateAngle(tb[10], tb[11], M_PI_2, 1);
    ta[1] = KP_CreateAngle(tb[11], tb[10], M_PI_2 * (-1), 1);
    KG_AddBallBindedTexture( 4, 1, tb[10], 2*M_PI_2, 13, tb[10], 0, 13, tb[11],0, 13, tb[11], 2*M_PI_2, 13);


    tb[9]->fx = -2;

    KP_SetGravity(0, 4);

    this_thread::sleep_for(chrono::milliseconds(500));

	while (1)
	{

		glClear(GL_COLOR_BUFFER_BIT);

        KP_DoGravity();
		KP_UpdatePhysics();

		// t4 = KA_TestRuntime(&t1, &t2, &t3);

		KG_SetTexture(1);
		KG_Drawtext(to_string(KP_GetDeltaTime()), 10, 10, 15, 30, 14);
		KG_Drawtext(to_string(KA_GetTicks()), 10, 30, 15, 30, 14);
		KG_Drawtext(to_string(t3), 10, 330, 15, 30, 14);
		//KG_Drawtext(to_string(KA_TestRuntime()), 10, 360, 15, 30, 14);
        KG_Drawtext(to_string(CLOCKS_PER_SEC), 10, 380, 15, 30, 14);
        //KG_Drawtext(to_string((double)KA_TestRuntime()/(double)CLOCKS_PER_SEC), 10, 400, 15, 30, 14);

		KG_DrawEverything();

		KG_SwapWindow();

		lastmx = mx;
		lastmy = my;
		SDL_GetMouseState(&mausx, &mausy);
		mx = (float)mausx;
		my = (float)mausy;

		while (SDL_PollEvent(&theEvent)){
            if (theEvent.type == SDL_QUIT) {
                exit(0);
            }
        }
	}

}
