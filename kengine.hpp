#ifndef _ke_graphics_hpp_
#define _ke_graphics_hpp_


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <list>
#include <math.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>

using namespace std;


struct KE_Material{
    double rest;        // restitution (0-1) Sprungkraft
    double fric;        // friction (Reibung)
    unsigned colType;        // If bitwise AND = True => Collision
    unsigned ignType;        // If bitwise AND = False => Collision
    unsigned props;          // True False properties:
                            // 1.Bit => True => Gravity Sensitiv
};

struct KE_Edge;
struct KE_Angle;

struct KE_Ball{
    double x;			// x of the basis (balance point)
    double y;			// y of the basis
    double r;			// rotation in radian
    double vx = 0;			// x-speed
    double vy = 0;			// y-speed
    double vr = 0;			// rotation-speed
    double fx = 0;			// x-force
    double fy = 0;			// y-force
    double fr = 0;			// rotation-force
    double ra;          // radius
    double m;		    // mass
    double c;           // circumcer.. (Umfang)
    KE_Material *mat;     // some atributes
    KE_Ball *colDone1 = nullptr;  // Ball collided with Edge => must not
    KE_Ball *colDone2 = nullptr;  //  collide with connected Balls anymore!
    list<KE_Edge*> edges;
    list<KE_Angle*> angles;
    KE_Ball *next;	    // next list element
};

struct KE_Edge{
    KE_Ball *ball1;
    KE_Ball *ball2;
    double len;           // Leangth
    double lenSqd;
    double massRate;       // Masse Ball1 im Verhaltnis zu Masse beider Balle
    double pow;           // Power between the balls (0-1)
    double siz;            // thickness
    KE_Material *mat;
    int solid = 1;      // can collide
    KE_Edge *next;
};

struct KE_Angle{
    KE_Ball *ball1;     // Angel at ball1
    KE_Ball *ball2;     // to ball 2
    double a;           // Angel
    double p;           // Power
    KE_Angle *next;
};


struct KE_BallBindSet{
    KE_Ball* ball[4];		// The 4 Edges are connected to Balls
    double angle[4];		// The angel manipulator |  with these you can displace the texture
    double distance[4];		// Distance manipulator	 |  edges from the original Ball center
    KE_BallBindSet *next;		// next List Element, 0 = End
};

struct KE_Image{
    int textureId;			// Id of the used texture
    int priority;			// We sort by this because near things have be drawn at last
    int type;               // 0 = with BallBindSet; 1 = like 0, but with angle, distance
    KE_BallBindSet *balls;
    KE_Image *next;
};

struct KE_BallSet{
    vector<KE_Ball*> ball;
    vector<KE_Edge*> edge;
    vector<KE_Angle*> angle;
};

struct KE_ColBall{
    KE_Ball *ball1;
    KE_Ball *ball2;
    KE_ColBall *next;
};

struct KE_ColEdge{
    int toHandle;       // 0 => do Ball Resolution; 1 => do Edge Resolution
    KE_Ball *ball;
    KE_Edge *edge;
    double rate;        // Collisionspunkt Entfernung zu ball1 im Verhältnis zur Distanz ball1 ball2
    double colx;
    double coly;
    KE_ColEdge *next;
};


// Functions:

// Allgemeine Funktionen:

void KA_InitAll(const string &title, int windowWidth, int windowHeight, int graphicsWidth, int graphicsHeight, int fullscreen);
SDL_Window* KA_GetWindow();

unsigned KA_GetTicks();
unsigned KA_TestRuntime(int *r1 = nullptr, int *r2 = nullptr, int *r3 = nullptr);


// Physik:

void KP_SetGravity(double gravityX, double gravityY);
void KP_DoGravity();

void KP_GetAllCollisions();
void KP_DoImpulseResolutions();
void KP_UpdatePhysics();
void KP_UpdatePositions();
void KP_UpdateEdges();
void KP_UpdateAngles();
void KP_UpdateDeltaTime();
void KP_ResetLastLoopTime();
void KP_WaitRemainingFrameTime();
double KP_GetDeltaTime();

void KP_NormalizeRelatedBalls(KE_Ball *ball);
void KP_NormalizeBall(KE_Ball *ball, KE_Edge *edge, KE_Ball *ball2);

KE_Ball* KP_CreateBall(double x, double y, double rotation, KE_Material *material, double radius, double m);

KE_Material* KP_CreateMaterial(double restitution, double friction, unsigned colType = 0xffffffff, unsigned ignoreType = 0, unsigned props = 0xffffffff);

KE_Edge* KP_CreateEdge(KE_Ball *ball1, KE_Ball *ball2, KE_Material *material, double power = 1);
KE_Edge* KP_CreateEdge(KE_Ball *ball1, KE_Ball *ball2, double power = 1);

KE_Angle* KP_CreateAngle(KE_Ball *ball1, KE_Ball *ball2, double angleValue, double power);

void KP_SetStandartMaterial(KE_Material *material);


// Grafik:

void KG_LoadTexturesWithTxdFile(const string &descriptionFile);
void KG_SetTexture(int textureId);
void KG_DrawSquare(GLfloat x1, GLfloat y1, GLfloat x2, GLfloat y2, GLfloat x3, GLfloat y3, GLfloat x4, GLfloat y4);
void KG_DrawBox(GLfloat x, GLfloat y, GLfloat w, GLfloat h);
void KG_SwapWindow();
void KG_ClearEverything();
void KG_Drawtext(const string &text, GLfloat x, GLfloat y, GLfloat width, GLfloat height, GLfloat distance);
void KG_DrawEverything();

GLubyte* KG_LoadPNGinArray(const string &file, unsigned* w, unsigned* h);

KE_Image* KG_AddBallBindedTexture(int textureId, int priority,
    KE_Ball* ball1, double angle1, double distance1, KE_Ball* ball2, double angle2, double distance2,
    KE_Ball* ball3, double angle3, double distance3, KE_Ball* ball4, double angle4, double distance4);
KE_Image* KG_AddBallBindedTexture(int textureId, int priority,
    KE_Ball* ball1, KE_Ball* ball2, KE_Ball* ball3, KE_Ball* ball4);
KE_Image* KG_AddBallBindedTexture(int textureId, int priority, KE_Ball* ball, double radius);



#endif
