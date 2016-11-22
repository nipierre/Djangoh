#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>

using namespace std;

//Structs

struct Events { vector<int> I;
                vector<string> particle;
                vector<int> KS;
                vector<int> KF;
                vector<int> orig;
                vector<double> p_x;
                vector<double> p_y;
                vector<double> p_z;
                vector<double> E;
                vector<double> m;
                double xBj;
                double y;
                double Q2; };


//Graphic Style

int fMarkerColor[2] = {4,2};
int fMarkerStyle[2] = {24,20};


//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);
