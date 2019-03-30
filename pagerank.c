#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include <string.h>

typedef struct CRS_data CRS_data;
typedef struct PageRank PageRank;
typedef struct TopScore TopScore;

struct CRS_data read_graph_from_file();
struct PageRank PageRank_Iterations();
struct TopScore top_n_webpages();

double dangling_nodes_W(int dangling_nodes[], double x_0[], int nodes);
double StoppingCriterion(int eps, double x_before[], double x_now[], int nodes);
double checkSign(double number);

void sortScore(double *x_k, int *node_ID, int start, int end, int nodes);
void mergeAll(double *x_k, int *node_ID, int start1, int end1, int start2, int end2, int nodes);

struct CRS_data
{
    _Bool error;
    int *arr_row_ptr;
    int *arr_col;
    double *arr_val;
    int *dangling_nodes;
    int trueEdges;
    int nodes;
};

struct PageRank
{
    double *x_k;
};

struct TopScore
{
    double *pagerankTopScore;
    int *node_number;
};

int main(void)
{
    struct CRS_data pageRankCRS;
    struct PageRank pageRankScore;
    struct TopScore highestScore;

    double eps, d;

    clock_t start, end;
    double seconds;
    int n;
    n = 100;
    d = 0.85;

    //d = 1;

    start = clock();

    pageRankCRS = read_graph_from_file(&pageRankCRS.arr_row_ptr, &pageRankCRS.arr_col, &pageRankCRS.arr_val, &pageRankCRS.dangling_nodes);
    end = clock();
    seconds = (double)(end - start) / CLOCKS_PER_SEC;

    printf("read from file took %f seconds \n", seconds);
    if (pageRankCRS.error)
    {
        return 1;
    }

    start = clock();
    pageRankScore = PageRank_Iterations(pageRankScore.x_k, pageRankCRS.arr_row_ptr, pageRankCRS.arr_col, pageRankCRS.arr_val, pageRankCRS.dangling_nodes, pageRankCRS.nodes, pageRankCRS.trueEdges, eps, d);
    end = clock();
    seconds = (double)(end - start) / CLOCKS_PER_SEC;

    printf("scoreiteration took %f seconds \n", seconds);

    start = clock();

    highestScore = top_n_webpages(pageRankScore.x_k, &highestScore.pagerankTopScore, &highestScore.node_number, pageRankCRS.nodes, n);

    end = clock();
    seconds = (double)(end - start) / CLOCKS_PER_SEC;

    printf("top_score took %f seconds \n", seconds);
    for (int l = 0; l < n; l++)
    {
        printf("l=%d score = %f with ", l, highestScore.pagerankTopScore[l]);
        printf("webgraph ID: %d \n", highestScore.node_number[l]);
    }

    free(pageRankCRS.arr_val);
    free(pageRankCRS.arr_col);
    free(pageRankCRS.arr_row_ptr);
    free(pageRankCRS.dangling_nodes);

    free(highestScore.pagerankTopScore);

    free(highestScore.node_number);

    //free(pageRankScore.x_k);
}

struct TopScore top_n_webpages(double *x_k, double *sorted_x, int *node_ID, int nodes, int n)
{

    printf("n = %d\n", n);
    sorted_x = malloc(n * sizeof *sorted_x);
    node_ID = malloc(n * sizeof *node_ID);
    int *x_k_node_ID = malloc(nodes * sizeof *x_k_node_ID);

    for (int i = 0; i < nodes; i++)
    {
        x_k_node_ID[i] = i;
    }

    sortScore(x_k, x_k_node_ID, 0, nodes, nodes);

    int teller = 0;
    int j = 0;
    for (int i = nodes; i > nodes - n; i--)
    {
        //printf("i = %d score %f \n", i, x_k[i]);
        sorted_x[j] = x_k[i];
        //printf("i = %d score %f \n", i, sorted_x[j]);
        node_ID[j] = x_k_node_ID[i];
        //printf("j = %d node %d \n", j, node_ID[j]);
        //teller++;
        //printf("teller %d \n", teller);
        j++;
    }

    TopScore score = {sorted_x, node_ID};

    return score;
}

void sortScore(double *x_k, int *node_ID, int start, int end, int nodes)
{
    // finner noden i midten:
    int middle;
    if (start < end)
    {
        middle = (end + start) / 2;
        // right_recursion
        sortScore(x_k, node_ID, start, middle, nodes);
        //left_recursion
        sortScore(x_k, node_ID, middle + 1, end, nodes);
        //merge all
        mergeAll(x_k, node_ID, start, middle, middle + 1, end, nodes);
    }
}

void mergeAll(double *x_k, int *node_ID, int start1, int end1, int start2, int end2, int nodes)
{
    double temp_x_k[nodes];
    double temp_node_ID[nodes];

    int i, j, k;

    i = start1;
    j = start2;
    k = 0;

    while (i <= end1 && j <= end2)
    {
        if (x_k[i] < x_k[j])
        {
            temp_x_k[k] = x_k[i];
            temp_node_ID[k] = node_ID[i];
            k++;
            i++;
        }
        else
        {
            temp_x_k[k] = x_k[j];
            temp_node_ID[k] = node_ID[j];

            k++;
            j++;
        }
    }

    while (i <= end1)
    {
        temp_x_k[k] = x_k[i];
        temp_node_ID[k] = node_ID[i];

        k++;
        i++;
    }

    while (j <= end2)
    {
        temp_x_k[k] = x_k[j];
        temp_node_ID[k] = node_ID[j];

        k++;
        j++;
    }

    int l = 0;
    for (i = start1; i <= end2; i++)
    {
        x_k[i] = temp_x_k[l];
        node_ID[i] = temp_node_ID[l];
        l++;
    }
}

struct CRS_data read_graph_from_file(int *arr_row_ptr, int *arr_col, double *arr_val, int *dangling_nodes)
{
    int trueEdges = 0;

    int row;
    int col;
    int teller = 0;

    int from, to;

    FILE *fp;

    //fp = fopen("test2.txt", "r");
    //fp = fopen("test100.txt", "r");

    fp = fopen("web-NotreDame.txt", "r");

    int n = 0; //100;
    char line[128];
    int nodes, edges;
    int lines = 0;

    if (fp == NULL)
    {
        CRS_data e = {1};
        puts("File could not open");
        return e;
    }

    fgets(line, sizeof line, fp);
    fgets(line, sizeof line, fp);

    fgets(line, sizeof line, fp);
    sscanf(line, "# Nodes: %d Edges: %d", &nodes, &edges);
    fgets(line, sizeof line, fp);

    arr_row_ptr = malloc((nodes + 1) * sizeof *arr_row_ptr);
    int *col_nr = malloc((nodes + 1) * sizeof *col_nr);
    int *row_nr = malloc((nodes) * sizeof *row_nr);

    dangling_nodes = malloc(nodes * sizeof *dangling_nodes);

    for (int i = 0; i < nodes; i++)
    {
        dangling_nodes[i] = 0;
    }

    arr_row_ptr[0] = 0;

    while (fgets(line, sizeof line, fp))
    {
        sscanf(line, "%d %d", &from, &to);

        if (from != to)
        {
            col_nr[from] += 1;
            arr_row_ptr[to + 1] = arr_row_ptr[to + 1] + 1;
            trueEdges++;
        }
    }
    fclose(fp);

    int sum = 0;
    for (int i = 0; i < nodes + 1; i++)
    {
        sum += arr_row_ptr[i];
        arr_row_ptr[i] = sum;
        //printf("sum = %d \n", sum);
    }

    int m;

    for (int i = 1; i < nodes + 1; i++)
    {
        m = arr_row_ptr[i] - arr_row_ptr[i - 1];
        if (m == 0)
        {
            dangling_nodes[i - 1] = 1;
            printf("%d \n", i);
        }
    }

    for (int i = 0; i < nodes; i++)
    {
        if (dangling_nodes[i] != 0)
        {
            printf("i=%d dangling_nodes = %d \n", i, dangling_nodes[i]);
        }
    }
    n = trueEdges;

    arr_val = malloc(n * sizeof *arr_val);
    arr_col = malloc(n * sizeof *arr_col);

    int *arr_row = malloc(n * sizeof *arr_row);

    for (int i = 0; i < nodes; i++)
    {
        row_nr[i] = 0;
    }

    fp = fopen("web-NotreDame.txt", "r");
    //fp = fopen("test100.txt", "r");

    //fp = fopen("test2.txt", "r");

    if (fp == NULL)
    {
        CRS_data e = {1};
        puts("File could not open");
        return e;
    }

    lines = 0;
    int index_CRS = 0;
    int to_teller = 0;
    int funnet;

    while (fgets(line, sizeof line, fp))
    {
        if (lines > 3)
        {
            sscanf(line, "%d %d", &from, &to);
            if (from != to)
            {
                index_CRS = arr_row_ptr[to] + row_nr[to];
                if (row_nr[to] == 0)
                {
                    arr_col[index_CRS] = from;
                    arr_val[index_CRS] = 1.0 / col_nr[from];
                }
                else if (from > arr_col[index_CRS - 1])
                {
                    arr_col[index_CRS] = from;
                    arr_val[index_CRS] = 1.0 / col_nr[from];
                }
                else
                {
                    for (int i = 1; i <= row_nr[to]; i++)
                    {
                        if (from > arr_col[index_CRS - i])
                        {
                            funnet = i;
                            break;
                        }
                    }
                    for (int i = 0; i < funnet; i++)
                    {
                        int leftElement = i + 1;
                        arr_col[index_CRS - i] = arr_col[index_CRS - leftElement];
                        arr_val[index_CRS - i] = arr_val[index_CRS - leftElement];
                    }
                    arr_col[index_CRS - funnet] = from;
                    arr_val[index_CRS - funnet] = 1.0 / col_nr[from];
                }
                row_nr[to] += 1;
            }
        }
        lines++;
    }
    fclose(fp);

    CRS_data crs = {false, arr_row_ptr, arr_col, arr_val, dangling_nodes, trueEdges, nodes};

    free(col_nr);
    free(row_nr);

    free(arr_row);

    return crs;
}

struct PageRank PageRank_Iterations(double *x_k, int *arr_row_ptr, int *arr_col, double *arr_val, int *dangling_nodes, int nodes, int trueEdges, double eps, double d)
{
    double *x_0 = malloc(nodes * sizeof *x_0);
    x_k = malloc(nodes * sizeof *x_k);

    // inialize startvalue 1/nodes
    double x_0_start_value = 1 / (double)nodes;
    double b = 1 - d;
    int link = 0;
    int r = 0;
    double z, w;

    #pragma omp parallel for
    for (int i = 0; i < nodes; i++)
    {
        x_0[i] = x_0_start_value;
    }

    //Finner dangling_node_vekt
    w = dangling_nodes_W(dangling_nodes, x_0, nodes);

    // Finner x_1 :
    for (int l = 1; l < (nodes + 1); l++)
    {
        r = arr_row_ptr[l] - arr_row_ptr[l - 1];
        z = 0;
        for (int j = link; j < (link + r); j++)
        {
            // Finner produktet for rett rad og kolonne nr.
            z += arr_val[j] * x_0[arr_col[j]];
        }

        x_k[l - 1] = ((x_0_start_value * (b + (d * w))) + d * z);
        link += r;
    }

    eps = 0.000001;

    _Bool iterate = 1;

    double stopping;

    int iterations = 0;

    while (iterate)
    {
        iterate = 0;
        stopping = StoppingCriterion(eps, x_0, x_k, nodes);

        if (stopping < eps)
        {
            printf("PageRank score found \n");
        }
        else
        {
            link = 0;
            iterate = 1;

            // setter x_k til x_0 og finner neste x_k (x_k+1)
            for (int l = 0; l < nodes; l++)
            {
                x_0[l] = x_k[l];
            }

            w = dangling_nodes_W(dangling_nodes, x_0, nodes);

            for (int p = 1; p < nodes + 1; p++)
            {
                r = arr_row_ptr[p] - arr_row_ptr[p - 1];
                z = 0;
                for (int j = link; j < link + r; j++)
                {
                    z += arr_val[j] * x_0[arr_col[j]]; 
                }
                x_k[p - 1] = (x_0_start_value * (b + (d * w))) + d * z;
                link += r;
            }
        }
    }
    iterations++;
    printf("x_0");
    free(x_0);
    printf("x_k");

    PageRank score = {x_k};

    return score;
}

double dangling_nodes_W(int dangling_nodes[], double x_0[], int nodes)
{
    double w = 0;
    int ii;

    for (ii = 0; ii < nodes; ii++)
    {
        //printf("dangling %d \n", dangling_nodes[i]);
        if (dangling_nodes[ii] != 0)
        {
            // Find W_(k-1)
            w += x_0[ii];
            //printf("w = %f \n", w);
        }
    }
    return w;
}

// Funksjon som returnerer diff mellom x_k - x_(k-1)
double StoppingCriterion(int eps, double x_before[], double x_now[], int nodes)
{
    double stop;
    for (int i = 0; i < nodes; i++)
    {
        stop = x_now[i] - x_before[i];
        stop = checkSign(stop);
        // stop større enn eps. hopp ut og fortsett å iterere ny pagerank.
        if (stop > eps)
        {
            return stop;
        }
    }
    // alle elementer er mindre enn eps som betyr at siste pagerank score er funnet.
    return stop;
}

double checkSign(double number)
{
    if (number < 0)
    {
        number = (-1) * number;
    }
    return number;
} 