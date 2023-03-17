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
double Q = 500;
double antFactor = 0.8;
double randomFactor = 1.0;
int maxIterations = 20000;
#define RAND_FACT_START 1.00
#define RAND_FACT_END 0.10
#define RAND_FACT_LAST_ITER 2000
#define NUM_CITIES 10
#define EPSILON 0.0001
int best_tour_order[NUM_CITIES];
double probabilities[NUM_CITIES];
double graph[NUM_CITIES][NUM_CITIES];
double trails[NUM_CITIES][NUM_CITIES];
int current_index = 0;
double best_tour_length = -1;

void hamiltonian_cycle_graph( int v,
                              int e,
                              int dir_flag,
                              char* out_file,
                              char* ham_file );

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
        double length = graph[trail[NUM_CITIES - 1]][trail[0]];
        for (int i = 0; i < NUM_CITIES-1; i++)
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
            ant.visitCity(-1, rand() % NUM_CITIES);
        }
    }
    current_index = 0;
}

bool calculateProbabilities(Ant ant)
{
    int i = ant.trail[current_index];
    double pheromone = 0.0;
    for (int l = 0; l < NUM_CITIES; l++)
    {
        if (!ant.isVisited(l) && graph[i][l] > 0)
        {
            pheromone += pow(trails[i][l], alpha) * pow(1.0 / graph[i][l], beta) + EPSILON;
        }
    }
    if(pheromone == 0) {
        return false;
    }
    for (int j = 0; j < NUM_CITIES; j++)
    {
        if (ant.isVisited(j) || graph[i][j] == 0)
        {
            probabilities[j] = 0.0;
        }
        else
        {
            double numerator = pow(trails[i][j], alpha) * pow(1.0 / graph[i][j], beta) + EPSILON;
            probabilities[j] = max(EPSILON, (numerator) / (pheromone));
        }
    }
    return true;
}

int selectNextCity(Ant ant)
{
    int t = rand() % NUM_CITIES;
    if (randDouble() < randomFactor)
    {
        if ((!ant.isVisited(t)) && graph[ant.trail[current_index]][t] > 0) // ! should be inside
        {
            return t;
        }
    }
    if(!calculateProbabilities(ant)) {
        return -1;
    }
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
        ants.erase(std::remove_if(ants.begin(), ants.end(), 
                       [](Ant& ant) {
                         int nextCity = selectNextCity(ant);
                         if(nextCity == -1) {
                            return true;
                         }
                         ant.visitCity(current_index, nextCity);
                         return false;
                         }), ants.end());
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
    if (best_tour_order[0] == -1 && ants.size() > 0)
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
    int is_directed = 0;
    hamiltonian_cycle_graph(NUM_CITIES,(NUM_CITIES*3),is_directed,"output-graph.txt","ham-path.txt");
}

int main()
{
    // srand((unsigned)time(NULL));

    generateRandomMatrix();
    int numberOfAnts = (int)(NUM_CITIES * antFactor);
    for(int i = 0; i < NUM_CITIES; ++i) {
        probabilities[i] = 0.5;
        best_tour_order[i] = -1;
    }

    double best[maxIterations];
    printf("\n        initializing...\n");

    // for(int i = 0; i < NUM_CITIES; ++i) {
    //     for(int j = 0; j < NUM_CITIES; ++j) {
    //         if(graph[i][j] == 0) {
    //             graph[i][j] = randDouble() * 500 + 1;
    //         }
    //         // printf("%.2f ", graph[i][j]);
    //     }
    //     // printf("\n");
    // }

    for(int i = 0; i < maxIterations; ++i) {
        ants.clear();
        for(int j = 0; j < numberOfAnts; ++j) {
            ants.push_back(Ant());
        }
        setupAnts(numberOfAnts);
        randomFactor = max(RAND_FACT_END, RAND_FACT_START - i / RAND_FACT_LAST_ITER);
        moveAnts();
        updateTrails();
        updateBest();
        best[i] = best_tour_length;
    }

    printf("\n Best tour:");
    for(int& i : best_tour_order) {
        printf("%d -> ", i);
    }

    ofstream MyFile("results.txt");
    for(int i = 0; i < maxIterations; ++i) {
        MyFile << best[i];
        if(i != maxIterations-1)
            MyFile << ",";
    }
    MyFile.close();

    printf("%d\n", best_tour_order[0]);
    printf("Best tour length: %f", best_tour_length);
}

//TODO: Fix visiting city twice


/* This function writes a simple graph with a Hamiltonian cycle to one
   file and a Hamiltonian cycle to another file. The graph will
   have max(e,v) edges. The graph can be directed or undirected. It is
   assumed that e <= v(v-1)/2 if the graph is undirected, and that
   e <= v(v-1) if the grpah is directed. (In this program,
   this assured because of the call to fix_imbalanced_graph.)

   To generate a random graph with a
   Hamiltonian cycle, we begin with a random permutation
   ham[0],...,ham[v-1] (v = number of vertices).  We then generate edges

   (ham[0],ham[1]),(ham[1],ham[2]),...,(ham[v-2],ham[v-1]),(ham[v-1],ham[0]),

   so that the graph has the Hamiltonian cycle

                     ham[0],...,ham[v-1],ham[0].

   Finally,  we  add random edges to produce the desired number
   of edges.
 */
void print_graph( int v,
                  int e,
                  char* out_file,
                  int* adj_matrix,
                  int dir_flag );
                
void init_array( int* a, int end );
int ran( int k );
void swap( int* a, int *b );
void permute( int* a, int n );

void hamiltonian_cycle_graph( int v,
                              int e,
                              int dir_flag,
                              char* out_file,
                              char* ham_file )
{
   int i, j, k, l, count, index, *adj_matrix, *ham;
   FILE *fptr;

   if ( ( adj_matrix = ( int * ) calloc( v * v, sizeof( int ) ) )
        == NULL ) {
      printf( "Not enough room for this size graph\n" );
      return;
   }

   if ( ( ham = ( int * ) calloc( v, sizeof( int ) ) ) == NULL ) {
      printf( "Not enough room for this size graph\n" );
      free( adj_matrix );
      return;
   }

   printf( "\n\tBeginning construction of graph.\n" );

   /*  Generate a random permutation in the array ham. */
   init_array( ham, v );
   permute( ham, v );

   if ( ( fptr = fopen( ham_file, "w" ) ) == NULL ) {
      printf( "\n\t\t Could not open file %s.\n", ham_file );
      free( adj_matrix );
      free( ham );
      return;
   }

   /* print Hamiltonian cycle and store required edges */
   for ( i = 0; i < v; i++ ) {
      fprintf( fptr, "%5d\n", 1 + ham[ i ] );
      k = ham[ i ];
      l = ham[ ( i + 1 ) % v ];
      if ( k > l && !dir_flag )
         swap( &k, &l );
      adj_matrix[ k * v + l ] = 1;
   }

   fprintf( fptr, "%5d\n", 1 + ham[ 0 ] );
   fclose( fptr );
   free( ham );

   for ( count = v; count < e; ) {
      if ( ( i = ran( v ) ) == ( j = ran( v ) ) )
         continue;
      if ( i > j && !dir_flag )
         swap( &i, &j );
      index = i * v + j;
      if ( !adj_matrix[ index ] ) {
         adj_matrix[ index ] = 1;
         count++;
      }
   }

   print_graph( v, count, out_file, adj_matrix, dir_flag );
}

/* set a[ i ] = i, for i = 0,...,end - 1 */
void init_array( int* a, int end )
{
   int i;

   for ( i = 0; i < end; i++ )
      *a++ = i;
}

/* randomly permute a[ 0 ],...,a[ n - 1 ] */
void permute( int* a, int n )
{
   int i;

   for ( i = 0; i < n - 1; i++ )
      swap( a + i + ran( n - i ), a + i );
}

void swap( int* a, int *b )
{
   int temp;

   temp = *a;
   *a = *b;
   *b = temp;
}

/* Return a random integer between 0 and k-1 inclusive. */
int ran( int k )
{
   return rand() % k;
}


void print_graph( int v,
                  int e,
                  char* out_file,
                  int* adj_matrix,
                  int dir_flag )
{
   int i, j, index;
   FILE *fp;

   if ( ( fp = fopen( out_file, "w" ) ) == NULL ) {
      printf( "Unable to open file %s for writing.\n", out_file );
      return;
   }
   printf( "\n\tWriting graph to file %s.\n", out_file );

   fprintf( fp, "%5d   %5d\n", v, e );

   if ( !dir_flag )
      for ( i = 1; i < v; i++ )
         for ( j = i + 1; j <= v; j++ ) {
            index = ( i - 1 ) * v + j - 1;
            if ( adj_matrix[ index ] ){
                double x = randDouble() * 500;
               fprintf( fp, "%5d   %5d   %5f\n", i, j, (adj_matrix[ index ] * (x)) );
               graph[i-1][j-1] = adj_matrix[ index ] * (x);
            }else{
                graph[i-1][j-1] = 0;
            }
         }
   else
      for ( i = 1; i <= v; i++ )
         for ( j = 1; j <= v; j++ ) {
            index = ( i - 1 ) * v + j - 1;
            if ( adj_matrix[ index ] ){
                double x = randDouble() * 500;
               fprintf( fp, "%5d   %5d   %5f\n", i, j, (adj_matrix[ index ] * (x)) );
               graph[i-1][j-1] = adj_matrix[ index ] * (x);
            }else{
                graph[i-1][j-1] = 0;
            }
         }
   fclose( fp );
   printf( "\tGraph is written to file %s.\n", out_file );
}