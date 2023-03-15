#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <time.h>

using namespace std;

double c = 1.0;
double alpha = 1;
double beta = 5;
double evaporation = 0.5;
double Q = 500;
double antFactor = 0.8;
double randomFactor = 0.01;
int maxIterations = 2000;
#define NUM_CITIES 20
#define TRAIL_SIZE 20
double probabilities[NUM_CITIES];
double graph[NUM_CITIES][NUM_CITIES];
double trails[NUM_CITIES][NUM_CITIES];
int current_index = 0;
int best_tour_order[NUM_CITIES];
int best_tour_length = -1;

class Ant
{
public:
    int visited[NUM_CITIES];
    int trail[TRAIL_SIZE];

    int trailSize() {
        for(int i = 0; i < TRAIL_SIZE; ++i) {
            if(trail[i] == -1) {
                return i+1;
            }
        }
        return TRAIL_SIZE;
    }

    void visitCity(int current_index, int city)
    {
        trail[current_index + 1] = city;
        visited[city] = true;
    }

    bool isVisited(int i)
    {
        return visited[i];
    }

    double trailLength()
    {
        double length = graph[trail[trailSize() - 1]][trail[0]];
        for (int i = 0; i < trailSize() - 1; i++)
        {
            length += graph[trail[i]][trail[i + 1]];
        }
        return length;
    }

    void clear()
    {
        for (int i = 0; i < TRAIL_SIZE; ++i)
        {
            trail[i] = -1;
        }
    }
};

vector<Ant> ants;

double randDouble() {
    return ((double)rand()/(double)RAND_MAX);
}

void setupAnts(int number_of_ants)
{
    for (int i = 0; i < number_of_ants; ++i)
    {
        for (Ant ant : ants)
        {
            ant.clear();
            ant.visitCity(-1, rand() % NUM_CITIES);
        }
    }
    current_index = 0;
}

int selectNextCity(Ant ant)
{
    int t = rand() % (NUM_CITIES - current_index);
    if (rand() < randomFactor)
    {
        if (!(ant.isVisited(t)))
        {
            return t;
        }
    }
    calculateProbabilities(ant);
    double r = randDouble();
    double total = 0;
    for (int i = 0; i < NUM_CITIES; i++) {
        total += probabilities[i];
        if (total >= r) {
            return i;
        }
    }
}

void moveAnts()
{
    for (int i = current_index; i < NUM_CITIES; ++i)
    {
        for (Ant ant : ants)
        {
            ant.visitCity(current_index, selectNextCity(ant));
        }
        current_index++;
    }
}

void calculateProbabilities(Ant ant)
{
    int i = ant.trail[current_index];
    double pheromone = 0.0;
    for (int l = 0; l < NUM_CITIES; l++)
    {
        if (!ant.isVisited(l))
        {
            pheromone += pow(trails[i][l], alpha) * pow(1.0 / graph[i][l], beta);
        }
    }
    for (int j = 0; j < NUM_CITIES; j++)
    {
        if (ant.isVisited(j))
        {
            probabilities[j] = 0.0;
        }
        else
        {
            double numerator = pow(trails[i][j], alpha) * pow(1.0 / graph[i][j], beta);
            probabilities[j] = numerator / pheromone;
        }
    }
}

void updateTrails()
{
    for (int i = 0; i < NUM_CITIES; i++)
    {
        for (int j = 0; j < NUM_CITIES; j++)
        {
            trails[i][j] *= evaporation;
        }
    }
    for (Ant a : ants)
    {
        double contribution = Q / a.trailLength();
        for (int i = 0; i < NUM_CITIES - 1; i++)
        {
            trails[a.trail[i]][a.trail[i + 1]] += contribution;
        }
        trails[a.trail[NUM_CITIES - 1]][a.trail[0]] += contribution;
    }
}

void clone(int n, int* src, int* dest) {
    for(int i = 0; i < n; ++i) {
        dest[i] = src[i];
    }
}

void updateBest()
{
    if (best_tour_order[0] == -1)
    {
        clone(TRAIL_SIZE, ants[0].trail, best_tour_order);
        best_tour_length = ants[0].trailLength();
    }
    for (Ant a : ants)
    {
        if (a.trailLength() < best_tour_length)
        {
            best_tour_length = a.trailLength();
            clone(TRAIL_SIZE, a.trail, best_tour_order);
        }
    }
}

void generateRandomMatrix() {
    for(int i = 0; i < NUM_CITIES; ++i) {
        for(int j = 0; j < NUM_CITIES; ++j) {
            graph[i][j] = randDouble();
        }
    }
}

int main()
{
    srand((unsigned)time(NULL));

    generateRandomMatrix();
    int numberOfAnts = (int)(NUM_CITIES * antFactor);
    setupAnts(numberOfAnts);
    for(int i = 0; i < numberOfAnts; ++i) {
        ants.push_back(Ant());
    }

    for(int i = 0; i < maxIterations; ++i) {
        moveAnts();
        updateTrails();
        updateBest();
    }
}