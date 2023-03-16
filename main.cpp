#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <algorithm>

using namespace std;

double c = 1.0;
double alpha = 1;
double beta = 5;
double evaporation = 0.5;
double Q = 50;
double antFactor = 0.8;
double randomFactor = 0.01;
int maxIterations = 1000;
#define RAND_FACT_START 1.00
#define RAND_FACT_END 0.05
#define RAND_FACT_LAST_ITER 500
#define NUM_CITIES 20
double probabilities[NUM_CITIES];
double graph[NUM_CITIES][NUM_CITIES];
double trails[NUM_CITIES][NUM_CITIES];
int current_index = 0;
int best_tour_order[NUM_CITIES];
double best_tour_length = -1;

class Ant
{
public:
    bool visited[NUM_CITIES];
    int trail[NUM_CITIES];

    Ant() {
    }

    int trailSize() {
        for(int i = 0; i < NUM_CITIES; ++i) {
            if(trail[i] == -1) {
                return i+1;
            }
        }
        return NUM_CITIES;
    }

    void visitCity(int current_index, int city)
    {
        trail[current_index+1] = city;
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
        for (int i = 0; i < NUM_CITIES; ++i)
        {
            trail[i] = -1;
            visited[i] = false;
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
        for (Ant& ant : ants)
        {
            ant.clear();
            ant.visitCity(-1, rand() % 20 + 1);
        }
    }
    current_index = 0;
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
            // printf("\n phero: %f", pow(trails[i][l], alpha));
        }
    }
    int c = 0;
    for (int j = 0; j < NUM_CITIES; j++)
    {
        if (ant.isVisited(j))
        {
            probabilities[j] = 0.0;
            c++;
        }
        else
        {
            double numerator = pow(trails[i][j], alpha) * pow(1.0 / graph[i][j], beta);
            probabilities[j] = (numerator+0.01) / (pheromone+0.01);
        }
    }
    // if(c == NUM_CITIES) {
        // printf("wtfau %d %d", current_index, ant.trailSize());
    // }
}

int selectNextCity(Ant ant)
{
    int t = rand() % (NUM_CITIES - current_index);
    if (randDouble() < randomFactor)
    {
        if (!(ant.isVisited(t)))
        {
            return t;
        }
    }
    calculateProbabilities(ant);
    double r = randDouble();
    double mx = 0;
    int mxidx = 0;
    for(int i = 0; i < NUM_CITIES; ++i) {
        if(probabilities[i] > mx) {
            mx = probabilities[i];
            mxidx = i;
        }
    }
    return mxidx;
}

void moveAnts()
{
    for (int i = current_index; i < NUM_CITIES-1; ++i)
    {
        for (Ant& a: ants)
        {
            a.visitCity(current_index, selectNextCity(a));
        }
        current_index++;
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
    for (Ant& a : ants)
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
        clone(NUM_CITIES, ants[0].trail, best_tour_order);
        best_tour_length = ants[0].trailLength();
    }
    for (Ant& a : ants)
    {
        if (a.trailLength() < best_tour_length)
        {
            best_tour_length = a.trailLength();
            clone(NUM_CITIES, a.trail, best_tour_order);
        }
    }
}

void generateRandomMatrix() {
    for(int i = 0; i < NUM_CITIES; ++i) {
        for(int j = 0; j < NUM_CITIES; ++j) {
            graph[i][j] = randDouble() * 500;
        }
    }
}

int main()
{
    srand((unsigned)time(NULL));

    generateRandomMatrix();
    int numberOfAnts = (int)(NUM_CITIES * antFactor);
    for(int i = 0; i < numberOfAnts; ++i) {
        ants.push_back(Ant());
    }
    for(int i = 0; i < NUM_CITIES; ++i) {
        probabilities[i] = 0.5;
        best_tour_order[i] = -1;
    }

    for(int i = 0; i < maxIterations; ++i) {
        setupAnts(numberOfAnts);
        randomFactor = max(RAND_FACT_END, RAND_FACT_START - i / RAND_FACT_LAST_ITER);
        moveAnts();
        updateTrails();
        updateBest();
    }

    printf("\n Best tour:");
    for(int& i : best_tour_order) {
        printf("%d -> ", i);
    }
    printf("%d\n", best_tour_order[0]);
    printf("Best tour length: %f", best_tour_length);
}