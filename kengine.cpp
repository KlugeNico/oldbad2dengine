//////////////////////////////////////////////////////////////////////
//////////////    Name:          Kluge Engine           //////////////
//////////////    Autor:         Nico Kluge             //////////////
//////////////    Begonnen:      Oktober 2014           //////////////
//////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
#include <thread>
#include <chrono>
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdint.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
#include "kengine.hpp"
#include "lodepng.h"

using namespace std;


GLuint *KE_textures;
SDL_Window* KE_window;
SDL_Event KE_event;

#define TEXTURES_MAX 10000



double KP_gravityX;
double KP_gravityY;

int KP_frameTime;

double dt;       // delta Time of last run for correct physic update
timeval KP_lastLoopTime;

KE_Ball *KE_firstBall;
KE_Edge *KE_firstEdge;
KE_Angle *KE_firstAngle;
KE_ColBall *KE_firstColBall;
KE_ColEdge *KE_firstColEdge;
KE_Image *KE_firstImage;

KE_Material *KE_standartMaterial;


// Functions:

void KA_InitAll(const string &title, int windowWidth, int windowHeight, int graphicsWidth, int graphicsHeight, int fullscreen){
	SDL_Init(SDL_INIT_EVERYTHING);
	atexit(SDL_Quit);
	if (!fullscreen){
		KE_window = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
			windowWidth, windowHeight, SDL_WINDOW_OPENGL);
	}
	else {
		KE_window = SDL_CreateWindow(title.c_str(), SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
			windowWidth, windowHeight, SDL_WINDOW_FULLSCREEN | SDL_WINDOW_OPENGL);
	}
	SDL_GLContext glcontext = SDL_GL_CreateContext(KE_window);
	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	KE_textures = new GLuint[TEXTURES_MAX];
	glGenTextures(TEXTURES_MAX, KE_textures);
	glClearColor(1.0, 0.0, 0.5, 1.0); // Farbe zum Löschen setzen
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, graphicsWidth, graphicsHeight, 0, -1, 1);
	KP_gravityX = 0;
    KP_gravityY = 0;
    dt = 0;
    KP_lastLoopTime.tv_sec = 0;
    KP_frameTime = 1000 / 50;
	KE_firstBall = nullptr;
	KE_firstEdge = nullptr;
	KE_firstAngle = nullptr;
    KE_firstColBall = nullptr;
    KE_firstColEdge = nullptr;
    KE_firstImage = nullptr;
}

GLubyte *KG_LoadPNGinArray(const string &file, unsigned* w, unsigned* h){
    vector<unsigned char> image;
    unsigned error = lodepng::decode(image, *w, *h, file.c_str());
	//If the loading went ok, convert to array
	if (!error){
		GLubyte *imagePtr;
		GLubyte *pixels;
		int arrayLeangth = *w * *h * 4;
		pixels = new GLubyte[arrayLeangth];	// Point new created array
		for (int i = 0; i < arrayLeangth; i++){
			*pixels = image[i];		// Copy the Pixelbuffer
			*(pixels++);			// Maybe I could use std::copy instead
		}
		pixels = pixels - arrayLeangth;
		return pixels;
	}
	else {
		cout << "Error: by reading PNG File: " << file.c_str() << endl;
		cout << "Error Text: " << lodepng_error_text(error) << endl;
		return nullptr;
	}
}

void KG_LoadTexturesWithTxdFile(const string &descriptionFile){
	ifstream descrFile;	// The Descriptionfile we work with
	descrFile.open(descriptionFile.c_str(), fstream::binary);
	string pngFile;		// Texture Filename
	int exit;			// Exit while(1) loop
	char curSign;		// Current Sign we readed from descrFile
	int id;
	unsigned w, h;
	GLubyte *pixels = nullptr;

	while (1){
		pngFile = "";
		exit = 0;
		curSign = '0';
		id = 0;

		while (curSign != '$'){		// Until a new Descriptionline comes
			if (descrFile.eof()){ exit = 1; break; }	// Finish if eof
			descrFile.read(&curSign, 1);	// Read a char
		}
		if (exit){ break; }

		// Read the ID from Descriptionfile
		for (int i = TEXTURES_MAX / 10; i > 0; i = i / 10){	// sum up char ASCII to int
			descrFile.read(&curSign, 1);
			id = id + (curSign & 0x00F) * i;	// Make char ASCII to int
		}
		descrFile.read(&curSign, 1);	// read not used char
		descrFile.read(&curSign, 1);	// and next not used char (normaly '"')
		descrFile.read(&curSign, 1);	// read first Letter
		while (curSign != '"'){			// until Filename ends ('"')
			pngFile.append(1, curSign);		// append Letter to String
			descrFile.read(&curSign, 1);
		}

		pixels = KG_LoadPNGinArray(pngFile, &w, &h);

		glBindTexture(GL_TEXTURE_2D, *(KE_textures + id));
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, 4, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

		printf("Loadtexture id %4d  w %4d  h %4d  File: \"", id, w, h);
		printf(pngFile.c_str());	// print everything to console
		printf("\"\n");
	}
	descrFile.close();
}



void KG_ClearEverything(){
	glClear(GL_COLOR_BUFFER_BIT);
}

void KG_SwapWindow(){
	SDL_GL_SwapWindow(KE_window);
}

SDL_Window* KA_GetWindow(){
	return KE_window;
}


KE_Image* KG_AddBallBindedTexture(int textureId, int priority,
		KE_Ball* ball1, double angle1, double distance1, KE_Ball* ball2, double angle2, double distance2,
		KE_Ball* ball3, double angle3, double distance3, KE_Ball* ball4, double angle4, double distance4){

	KE_Image *texture = new KE_Image;
	KE_BallBindSet *ballBindSet = new KE_BallBindSet;
	texture->balls = ballBindSet;
	texture->priority = priority;

    // We have to sort them by priority, because we want to draw them in the right order
	if (KE_firstImage){
		KE_Image *texture2 = KE_firstImage;
		if (texture->priority <= texture2->priority){
			texture->next = texture2;
			KE_firstImage = texture;
		}
		else {
			int done = 0;
			while (texture2->next){
				if (texture->priority <= texture2->next->priority){
					texture->next = texture2->next;
					texture2->next = texture;
					done = 1;
					break;
				}
				else {
					texture2 = texture2->next;
				}
			}
			if (!done){
				texture2->next = texture;
				texture->next = nullptr;
			}
		}
	}
	else {
		KE_firstImage = texture;
		texture->next = nullptr;
	}

    texture->type = 0;
	texture->textureId = textureId;
	texture->balls->ball[0] = ball1;
    texture->balls->ball[1] = ball2;
	texture->balls->ball[2] = ball3;
	texture->balls->ball[3] = ball4;

    if(angle1||distance1||angle2||distance2||angle3||distance3||angle4||distance4){
        texture->type = 1;
        texture->balls->angle[0] = angle1;
        texture->balls->angle[1] = angle2;
        texture->balls->angle[2] = angle3;
        texture->balls->angle[3] = angle4;
        texture->balls->distance[0] = distance1;
        texture->balls->distance[1] = distance2;
        texture->balls->distance[2] = distance3;
        texture->balls->distance[3] = distance4;
    }

	return texture;
}


KE_Image* KG_AddBallBindedTexture(int textureId, int priority,
		KE_Ball* ball1, KE_Ball* ball2,	KE_Ball* ball3, KE_Ball* ball4){
    return KG_AddBallBindedTexture(textureId, priority, ball1, 0, 0, ball2, 0, 0, ball3, 0, 0, ball4, 0, 0);
}

KE_Image* KG_AddBallBindedTexture(int textureId, int priority, KE_Ball* ball, double radius){
    radius = sqrt(radius*radius+radius*radius);
    return KG_AddBallBindedTexture(textureId, priority, ball, 1.25*M_PI, radius, ball, 1.75*M_PI, radius, ball, 0.25*M_PI, radius, ball, 0.75*M_PI, radius);
}


void KG_SetTexture(int textureId){
	glBindTexture(GL_TEXTURE_2D, *(KE_textures + textureId));
}


void KG_DrawSquare(GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2, GLfloat x3, GLfloat y3, GLfloat x4, GLfloat y4){
	glBegin(GL_QUADS);
	glTexCoord2f(0, 0);
	glVertex2f(x1, y1);

	glTexCoord2f(1, 0);
	glVertex2f(x2, y2);

	glTexCoord2f(1, 1);
	glVertex2f(x3, y3);

	glTexCoord2f(0, 1);
	glVertex2f(x4, y4);
	glEnd();
}

void KG_DrawBox(GLfloat x, GLfloat y, GLfloat w, GLfloat h){
	glBegin(GL_QUADS);
	glTexCoord2f(0, 0);
	glVertex2f(x, y);

	glTexCoord2f(1, 0);
	glVertex2f(x+w, y);

	glTexCoord2f(1, 1);
	glVertex2f(x+w, y+w);

	glTexCoord2f(0, 1);
	glVertex2f(x, y+w);
	glEnd();
}

// Before this fuction we have to load a Texture with letters.
// 16 Signs per row and 8 Signs per column (in ASCII order)
void KG_Drawtext(const string &text, GLfloat x, GLfloat y, GLfloat width, GLfloat height, GLfloat distance){
	const char *ctext = text.c_str();

	GLfloat lw = 0.0625;		// Letter Width
	GLfloat lh = 0.125;		// Letter Height
	GLfloat tx = 0;			// Texture x start of the Sign
	GLfloat ty = 0;			// Texture y start of the Sign
	int signx = 0;
	int signy = 0;

	glBegin(GL_QUADS);

	for (int i = 0; i < text.length(); i++){
		signx = (int)ctext[i] % 16;
		signy = (int)ctext[i] / 16;
		tx = (float)signx * lw;
		ty = (float)signy * lh;

		glTexCoord2f(tx, ty);
		glVertex2f(x, y);

		glTexCoord2f(tx + lw, ty);
		glVertex2f(x + width, y);

		glTexCoord2f(tx + lw, ty + lh);
		glVertex2f(x + width, y + height);

		glTexCoord2f(tx, ty + lh);
		glVertex2f(x, y + height);

		x = x + distance;
	}

	glEnd();
}


void KG_DrawEverything(){
	int i = 0;
	KE_Image *texture = KE_firstImage;
	while (texture){
		KG_SetTexture(texture->textureId);
		switch (texture->type){
		case 0:
			glBegin(GL_QUADS);
			for (i = 0; i < 4; i++){
				glTexCoord2f(((i + 1) & 0x2) >> 1, (i & 0x2) >> 1);	// 0(0,0) 1(1,0) 2(1,1) 3(0,1)
				glVertex2f(texture->balls->ball[i]->x, texture->balls->ball[i]->y);
			}
			glEnd();
			break;
		case 1:
			glBegin(GL_QUADS);
			for (i = 0; i < 4; i++){
				glTexCoord2f(((i + 1) & 0x2) >> 1, (i & 0x2) >> 1);	// 0(0,0) 1(1,0) 2(1,1) 3(0,1)
				glVertex2f(texture->balls->ball[i]->x + texture->balls->distance[i] * cos(texture->balls->angle[i] + texture->balls->ball[i]->r),
					texture->balls->ball[i]->y + texture->balls->distance[i] * sin(texture->balls->angle[i] + texture->balls->ball[i]->r));
			}
			glEnd();
			break;
		default: break;
		}
		texture = texture->next;
	}
}


KE_Ball* KP_CreateBall(double x, double y, double rotation, KE_Material *material, double radius, double m){
	KE_Ball *ball = new KE_Ball;

    ball->next = KE_firstBall;
    KE_firstBall = ball;

	ball->x = x;
	ball->y = y;
	ball->r = rotation;
    ball->mat = material;
	ball->ra = radius;
	ball->m = m;
    ball->c = 2*radius*M_PI;

	return ball;
}


KE_Edge* KP_CreateEdge(KE_Ball *ball1, KE_Ball *ball2, KE_Material *material, double power){
    KE_Edge *edge = new KE_Edge;

    edge->next = KE_firstEdge;
    KE_firstEdge = edge;

    edge->ball1 = ball1;
    edge->ball2 = ball2;
    edge->lenSqd = (ball1->x - ball2->x) * (ball1->x - ball2->x) + (ball1->y - ball2->y) * (ball1->y - ball2->y);
    edge->len = sqrt(edge->lenSqd);
    edge->siz = ball1->ra;
    edge->massRate = ball1->m / (ball1->m + ball2->m);
    edge->pow = power;
    edge->solid = 1;
    edge->mat = material;

    ball1->edges.push_front(edge);
    ball2->edges.push_front(edge);

    return edge;
}

KE_Edge* KP_CreateEdge(KE_Ball *ball1, KE_Ball *ball2, double power){
    return KP_CreateEdge(ball1, ball2, KE_standartMaterial, power);
}


KE_Angle* KP_CreateAngle(KE_Ball *ball1, KE_Ball *ball2, double angleValue, double power){
    KE_Angle *angle = new KE_Angle;

    angle->next = KE_firstAngle;
    KE_firstAngle = angle;

    angle->a = angleValue;
    angle->p = power;
    angle->ball1 = ball1;
    angle->ball2 = ball2;

    ball1->angles.push_front(angle);
    ball2->angles.push_front(angle);

    return angle;
}


KE_Material* KP_CreateMaterial(double restitution, double friction, unsigned colType, unsigned ignoreType, unsigned props){
	KE_Material *material = new KE_Material;

	material->rest = restitution;
	material->fric = friction;
	material->colType = colType;
	material->ignType = ignoreType;
	material->props = props;

	return material;
}

void KP_SetStandartMaterial(KE_Material *material){
    KE_standartMaterial = material;
}


void KP_NormalizeRelatedBalls(KE_Ball *ball){

    KE_Ball* workball;
    KE_Ball* editball;

    list<KE_Ball*> finishedBalls;
    list<KE_Ball*> toDoBalls;

    toDoBalls.push_front(ball);

    while(!toDoBalls.empty()){
        workball = toDoBalls.front();
        toDoBalls.pop_front();

        for (list<KE_Edge*>::iterator it = workball->edges.begin(); it != workball->edges.end(); it++){
            if((*it)->ball1 == workball){
                editball = (*it)->ball2;
            } else {
                editball = (*it)->ball1;
            }

            //editball-> TO DO

            toDoBalls.push_back(editball);
        }

    }

}


void KP_NormalizeBall(KE_Ball *ball, KE_Edge *edge, KE_Ball *ball2){

    //ball->x = ball2->x + TO DO

}


void KP_GetAllCollisions(){
	KE_Ball *ball = KE_firstBall;
    KE_Ball *ball1;
	KE_Ball *ball2;
	KE_Edge *edge = KE_firstEdge;
	KE_ColBall *colBall;
	KE_ColEdge *colEdge;
	double dSqd;		// Distance * Distance
    double disx;
    double disy;
    double dis1x;       // distance ball 1
    double dis1y;
    double dis2x;       // distance ball 2
    double dis2y;
    double dis1Sqd;
    double dis2Sqd;
    double deb;         // distance edge balls
    double deb1x;
    double deb1y;
    double deb2x;
    double deb2y;
    double debSqd;
    double ze;          // Zwischenergebnis
    double p;           // Proportion Distance betwenn col-pos and ball1 to deb
    double ara;         // avarage radius

    while (edge){
        while (ball){
            if (edge->ball1 != ball && edge->ball2 != ball){
                ball1 = edge->ball1;
                ball2 = edge->ball2;
                dis1x = ball->x - ball1->x;
                dis1y = ball->y - ball1->y;
                deb1x = ball2->x - ball1->x;
                deb1y = ball2->y - ball1->y;
                if ((dis1x*deb1x + dis1y*deb1y) > 0){
                    dis2x = ball->x - ball2->x;
                    dis2y = ball->y - ball2->y;
                    deb2x = deb1x * (-1);
                    deb2y = deb1y * (-1);
                    if ((dis2x*deb2x + dis2y*deb2y) > 0){
                        debSqd = deb1x*deb1x + deb1y*deb1y;
                        dis1Sqd = dis1x*dis1x + dis1y*dis1y;
                        dis2Sqd = dis2x*dis2x + dis2y*dis2y;
                        ze = debSqd + dis1Sqd - dis2Sqd;
                        dSqd = dis1Sqd - (ze*ze) / (4 * debSqd);
                        if (dSqd < (edge->siz+ball->ra) * (edge->siz+ball->ra)){
                            colEdge = new KE_ColEdge;
                            colEdge->next = KE_firstColEdge;
                            KE_firstColEdge = colEdge;
                            colEdge->ball = ball;
                            colEdge->edge = edge;
                            colEdge->rate = sqrt(dis1Sqd - dSqd) / sqrt(debSqd);
                            colEdge->colx = ball1->x + (ball2->x - ball1->x)*colEdge->rate;       // Kollisionspunkt
                            colEdge->coly = ball1->y + (ball2->y - ball1->y)*colEdge->rate;
                            ball->colDone1 = ball1;
                            ball->colDone2 = ball2;
                        }
                    }
                }
            }
        	ball = ball->next;
        }
        ball = KE_firstBall;
        edge = edge->next;
    }

    ball = KE_firstBall;
    ball2 = KE_firstBall;
	while (ball){
		while (ball2){
			if (ball != ball2){
                disx = ball2->x - ball->x;
                disy = ball2->y - ball->y;
				dSqd = disx * disx + disy * disy;
				if (dSqd < (ball->ra + ball2->ra) * (ball->ra + ball2->ra)){
                    colBall = new KE_ColBall;
                    colBall->next = KE_firstColBall;
                    KE_firstColBall = colBall;
                    colBall->ball1 = ball;
                    colBall->ball2 = ball2;
				}
			}
			ball2 = ball2->next;
		}
		ball2 = KE_firstBall;
		ball = ball->next;
	}
}

void KP_DoImpulseResolutions() {
	KE_ColBall *colBall = KE_firstColBall;
	KE_ColBall *nextColBall = nullptr;

    KE_ColEdge *colEdge = KE_firstColEdge;
    KE_ColEdge *nextColEdge = nullptr;

	double rvx;     // relative Geschwindigkeit nach x
	double rvy;     // relative Geschwindigkeit nach y
	double dx;     // Distanz x
    double dy;     // Distanz y
    double fx;      // Kraftvektor x
    double fy;      // Kraftvektor y
    double ze;      // Zwischenergebnis
    double dis;     // distance
    double disSqd;  // distance * distance
    double om;      // Mass Relation (Masseverhaeltnis)
    double om2;
    double zx;     // zwischenergebnis x
    double zy;     // zwischenergebnis y
    double friction; // Reibung
    double rr;      // Both Radius added
    double cor;     // Korrektur
    double colx;    // Kollisionspunkt
    double coly;
    double rate;    // Relation between Distance Ball1 and Ball to Ball1 Ball2
    double rest;    // Restitution (Sprungkraft)

    int i = 0;
    KE_Edge* edge;
    KE_Ball* ball;
    KE_Ball* ball1;
    KE_Ball* ball2;

    KG_SetTexture(1);

	while (colBall){
		ball = colBall->ball1;
        ball2 = colBall->ball2;

        // Ein paar Vorrechnungen:
        rvx = ball2->vx - ball->vx;     // Relative -v von ball zu ball2
        rvy = ball2->vy - ball->vy;
        dx = ball->x - ball2->x;       // Die Balldistanz
        dy = ball->y - ball2->y;
        disSqd = dx*dx+dy*dy;
        rest = (ball->mat->rest + ball2->mat->rest) / 2;

        om = ball2->m / (ball->m + ball2->m);   // Masseverhaeltnis
        rr = ball->ra + ball2->ra;

        // Der Ball wird abgestoßen:
        if ((dx*rvx + dy*rvy) > 0){
            ze = 2 * (dx*rvx+dy*rvy) / disSqd;    // maximale Laenge des Kraftvektors
            fx = dx*ze;                    // Kraftvektor auf max. Laenge
            fy = dy*ze;
            fx = fx*om*rest;       // Vektor abhängig von Masse und Sprungkraft
            fy = fy*om*rest;
            ball->fx = ball->fx + fx;
            ball->fy = ball->fy + fy;
        } else {
            // Bälle distanzieren:
            ze = ((rr*rr - disSqd) / disSqd) / 3;
            if(ze > 0.984375){ze = 0.984375;}
            ball->x = ball->x + dx*ze*om;
            ball->y = ball->y + dy*ze*om;
        }

        // Reibung:
        friction = ball->mat->fric + ball2->mat->fric;
        zx = dy; zy = -dx; cor = -1;
        if ((zx*rvx + zy*rvy) > 0){
            zx = -zx; zy = -zy; cor = 1;
        }
        ze = (zx*rvx+zy*rvy) / dis;
        zx = zx*ze;
        zy = zy*ze;
        ze = sqrt(zx*zx + zy*zy);
        if ((zx*zx + zy*zy) > friction*friction){
            zx = (zx * friction) / ze;
            zy = (zy * friction) / ze;
        }
        zx = zx - (zx/ze) * ball->vr * ball->c / (2*M_PI);
        zy = zy - (zy/ze) * ball->vr * ball->c / (2*M_PI);
        //ball->fx = ball->fx - zx;
        //ball->fy = ball->fy - zy;
        //ball->vr = cor * ze / ball->c;

		nextColBall = colBall->next;
		delete colBall;
		colBall = nextColBall;
	}


    while (colEdge){

        ball = colEdge->ball;
		ball1 = colEdge->edge->ball1;
		ball2 = colEdge->edge->ball2;
		edge = colEdge->edge;
        rate = colEdge->rate;
        colx = colEdge->colx;
        coly = colEdge->coly;

        rvx = (ball1->vx * rate + ball2->vx * (1-rate)) - ball->vx;     // Relative -v von ball zu edge
        rvy = (ball1->vy * rate + ball2->vy * (1-rate)) - ball->vy;

        dx = ball->x - colx;       // Die Balldistanz
        dy = ball->y - coly;
        dis = dx*dx+dy*dy;

        om = ball1->m*rate + ball2->m*(1-rate);
        om = om / (ball->m + om);       // Vektor abhängig von Masse
        rr = ball->ra + edge->siz;      // Soll Enfernung (Radius ball + Dicke Edge)

        rest = (ball->mat->rest + edge->mat->rest) / 2;      // Durchschnittliche Reibung und Sprungkraft
        friction = (ball->mat->fric + edge->mat->fric) / 2;

        // Alle Baelle werden abgestoßen:
        if ((dx*rvx + dy*rvy) > 0){
            ze = 2 * (dx*rvx+dy*rvy) / dis;    // maximale Laenge des Kraftvektors
            fx = dx*ze*rest;                         // Kraftvektor auf max. Laenge
            fy = dy*ze*rest;
            ball->fx = ball->fx + fx*om;
            ball->fy = ball->fy + fy*om;

            om2 = om - 1;
            ball1->fx = om2 * (1-rate) * fx;
            ball1->fy = om2 * (1-rate) * fy;
            ball2->fx = om2 * rate * fx;
            ball2->fy = om2 * rate * fy;
        }

        // Bälle distanzieren:
        ze = ((rr*rr - dis) / dis) / 3;
        if(ze > 0.984375){ze = 0.984375;}
        ball->x = ball->x + dx*ze*om;
        ball->y = ball->y + dy*ze*om;

        // Friction (Reibung)
        zx = dy; zy = -dx; cor = -1;
        if ((zx*rvx + zy*rvy) > 0){
            zx = -zx; zy = -zy; cor = 1;
        }
        ze = (zx*rvx+zy*rvy) / dis;
        zx = zx*ze;
        zy = zy*ze;
        ze = sqrt(zx*zx + zy*zy);
        //if ((zx*zx + zy*zy) > friction*friction){
        //    zx = (zx * friction) / ze;
        //    zy = (zy * friction) / ze;
        //}
        //zx = zx - (zx/ze) * ball->vr * ball->c / (2*M_PI);
        //zy = zy - (zy/ze) * ball->vr * ball->c / (2*M_PI);
        ball->fx = ball->fx + zx;
        ball->fy = ball->fy + zy;
        //ball->fr = ball->fr + cor * ze / ball->c;


        KG_Drawtext("EdgeCol "+to_string(rate), 10, 90 + i * 40, 10, 20, 10);

        nextColEdge = colEdge->next;
		delete colEdge;
		colEdge = nextColEdge;
        i++;
    }

	KE_firstColBall = nullptr;
    KE_firstColEdge = nullptr;
}


void KP_UpdateEdges(){

    KE_Edge* edge = KE_firstEdge;

    double dis;     // Real Distance between Balls
    double disx;
    double disy;
    double changex;     // Difference between Distance and Edgeleangth
    double changey;
    double rate;     // rate between real Distance and Edge Leangth

    while(edge){
        disx = edge->ball1->x - edge->ball2->x;
        disy = edge->ball1->y - edge->ball2->y;
        dis = sqrt(disx * disx + disy * disy);
        rate = edge->len / dis;
        changex = (disx * rate - disx) * edge->pow;
        changey = (disy * rate - disy) * edge->pow;

        edge->ball1->fx = edge->ball1->fx + changex * (1-edge->massRate);
        edge->ball1->fy = edge->ball1->fy + changey * (1-edge->massRate);
        edge->ball2->fx = edge->ball2->fx - changex * edge->massRate;
        edge->ball2->fy = edge->ball2->fy - changey * edge->massRate;

        edge = edge->next;
    }

}

void KP_UpdateAngles(){

    KE_Angle* angle = KE_firstAngle;

    double disx;        // Real Distance between Balls
    double disy;
    double rate;     // Difference between real Angle and Ballrotation => later change

    while(angle){
        disx = angle->ball1->x - angle->ball2->x;
        disy = angle->ball1->y - angle->ball2->y;
        rate = angle->ball1->r - atan2(disy, disx);
        if(rate > M_PI){rate = rate - 2 * M_PI;}
        if(rate < (-1 * M_PI)){rate = rate + 2 * M_PI;}
        rate = angle->ball1->fr + (angle->a - rate) * angle->p;
        angle->ball1->fr = rate;


        angle = angle->next;
    }

}

void KP_UpdatePositions(){
	KE_Ball *ball = KE_firstBall;
	while (ball){
        ball->vx = (ball->vx + ball->fx) * 0.995;
        ball->vy = (ball->vy + ball->fy) * 0.995;
        ball->vr = (ball->vr + ball->fr) * 0.95;
        ball->x = ball->x + ball->vx * dt;
        ball->y = ball->y + ball->vy * dt;
        ball->r = ball->r + ball->vr * dt;
        ball->fx = 0;
        ball->fy = 0;
        ball->fr = 0;
        if(ball->r > M_PI){ball->r = ball->r - 2 * M_PI;}
        if(ball->r < (-1 * M_PI)){ball->r = ball->r + 2 * M_PI;}
        ball = ball->next;
	}
}

void KP_SetGravity(double gravityX, double gravityY){
    KP_gravityX = gravityX;
    KP_gravityY = gravityY;
}

void KP_DoGravity(){
	KE_Ball *ball = KE_firstBall;
	while (ball){
        if (ball->mat->props & 0x00000001){
            ball->fx = ball->fx + KP_gravityX;
            ball->fy = ball->fy + KP_gravityY;
        }
        ball = ball->next;
	}
}

void KP_UpdateDeltaTime(){
    timeval time;
    unsigned div;
    gettimeofday(&time, 0);
    if(KP_lastLoopTime.tv_sec == 0) {
        div = 0;
    } else {
        div = (time.tv_sec - KP_lastLoopTime.tv_sec) * 1000000 + (time.tv_usec - KP_lastLoopTime.tv_usec);
    }
    KP_lastLoopTime = time;
    dt = double(div) / 1000000;
}

double KP_GetDeltaTime(){
    return dt;
}

void KP_ResetLastLoopTime(){
    KP_lastLoopTime.tv_sec = 0;
}

void KP_WaitRemainingFrameTime(){
    timeval time;
    int div;       // in Millisecs
    gettimeofday(&time, 0);
    div = (time.tv_sec - KP_lastLoopTime.tv_sec) * 1000 + (time.tv_usec - KP_lastLoopTime.tv_usec) / 1000;
    div = KP_frameTime - div;
    chrono::milliseconds dur(div);
    if(div > 0){
        this_thread::sleep_for(dur);
    }
}



void KP_UpdatePhysics(){
    KP_WaitRemainingFrameTime();
    KP_UpdateDeltaTime();
    KP_GetAllCollisions();
    KP_DoImpulseResolutions();
    KP_UpdateEdges();
    KP_UpdateAngles();
    KP_UpdatePositions();
}


// Return the Ticks since last call
unsigned KA_GetTicks(){
    unsigned timeSinceLastCall;
	static unsigned lastCall = 0;
    timeSinceLastCall = clock() - lastCall;
    lastCall = clock();
    return timeSinceLastCall;
}

unsigned KA_TestRuntime(int *r1, int *r2, int *r3){
	static int t1 = 0;
	static int t2 = 0;
	static int t3 = 0;
	static int runs = 0;
	static int firstTime = 0;

	double xDistance;
	double yDistance;
	double sqrdDistance;
	KE_Ball *ball = KE_firstBall;
	double x1 = ball->x;
	double y1 = ball->y;
	double x2 = 1.345363;
	double y2 = 1.14557;

	int sinceLastCall = clock() - firstTime;
	runs++;
	int a;
	firstTime = clock();
	for (int i = 0; i < 10000; i++){
		//FILL WITH Runtimetest 1:

		//END
	}
	t1 = t1 + clock() - firstTime;
	firstTime = clock();
	for (int i = 0; i < 10000; i++){
		//FILL WITH Runtimetest 2:
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		ball->x = ball->x;
		//END
	}
	t2 = t2 + clock() - firstTime;
	firstTime = clock();
	for (int i = 0; i < 10000; i++){
		//FILL WITH Runtimetest 3:
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		x1 = x1;
		//END
	}
	t3 = t3 + clock() - firstTime;

	if(r1){*r1 = t1 / runs;}
	if(r2){*r2 = t2 / runs;}
	if(r3){*r3 = t3 / runs;}

	firstTime = clock();
	return (unsigned)sinceLastCall;
}
