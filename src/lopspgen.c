/*
 This file was modified from the file plantri.c from the program plantri by Brendan McKay and Gunnar Brinkmann

 This is the copyright statement for plantri and associated utilities.

Copyright is jointly held by the authors
 Gunnar Brinkmann, University of Gent, gunnar.brinkmann@ugent.be
 Brendan McKay, Australian National University, brendan.mcKay@anu.edu.au

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this software except in compliance with the License.
A copy of the License is included in the package and you can also
view it at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

/*
Copyright 2023 Heidi Van den Camp

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
A copy of the License is included in the repository and you can also
view it at
        http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <stdio.h>
#include <stdlib.h>
#include<stdint.h>
#include <getopt.h>

#include <sys/times.h>
#if !defined(CLK_TCK) && !defined(_SC_CLK_TCK)
#include <unistd.h>
#endif
#if !defined(CLK_TCK) && defined(_SC_CLK_TCK)
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

#if !defined(CLK_TCK)
#define CLK_TCK 100
 /* If the CPU time stated by the program appears to be
    out by a constant ratio, the most likely explanation is
    that the code got to this point but 100 is the wrong guess. */
#endif

#define SWITCHES ":uolhm:c:sa"

static int       /* presence of command-line switches */
        oswitch,
        lswitch,
        sswitch,
        aswitch,
        uswitch;


#if CPUTIME
#include <sys/times.h>
#include <time.h>
#if !defined(CLK_TCK) && !defined(_SC_CLK_TCK)
#include <unistd.h>
#endif
#if !defined(CLK_TCK) && defined(_SC_CLK_TCK)
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

#if !defined(CLK_TCK)
#define CLK_TCK 100
 /* If the CPU time stated by the program appears to be
    out by a constant ratio, the most likely explanation is
    that the code got to this point but 100 is the wrong guess. */
#endif
#endif

#ifndef MAXN
#define MAXN 64            /* the maximum number of vertices;*/
#endif
#define MAXE (6*MAXN-12)   /* the maximum number of oriented edges */
#define MAXFACTOR 64        /* the maximum allowed inflation factor*/


typedef unsigned long long bigint;
static bigint nout[MAXFACTOR][6];      /* counts of output graphs, per connectivity */
static bigint nout_op[MAXFACTOR][6];   /* counts of output graphs, per connectivity, OP */

static bigint totalout;       /* Sum of nout[] (always) */
static bigint totalout_op;    /* Sum of nout_op[] (only if -o) */
static bigint totalquads;
static bigint total_quads_generated = 0;

typedef struct e /* The data type used for edges */
{
    int start;         /* vertex where the edge starts */
    int end;           /* vertex where the edge ends */
    int rightface;     /* face on the right side of the edge
                          note: only valid if make_dual() called */
    struct e *prev;    /* previous edge in clockwise direction */
    struct e *next;    /* next edge in clockwise direction */
    struct e *invers;  /* the edge that is inverse to this one */
    struct e *min;     /* the least of e and e->invers */
    int mark,index;    /* three ints for temporary use;
                          rf is only for the printing routines;
                          Only access mark via the MARK macros. */
} EDGE;

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE  1

static char *outfilename;  /* name of output file (NULL for stdout) */
static FILE *outfile;      /* output file for graphs */
static FILE *msgfile;      /* file for informational messages */

static int maxnv;          /* order of output graphs */
static int factor;         /* The desired inflation factor */
static int res,mod;        /* res/mod from command line (default 0/1) */
static int splitlevel,
        splitcount;     /* used for res/mod splitting */
static int min_connectivity;

static int needgroup;      /* Is group needed at end of scansimple()
                              and similar routines? */

static int gotone_nbop;    /* Used only by got_one() */
static int gotone_nbtot;

/* The variables below are used at each level of the iteration,
   updating and restoring as we move up and down the search tree */

static int nv;             /* number of vertices; they are 0..nv-1 */
static int ne;             /* number of directed edges (at most 6*nv-12) */

#define NUMEDGES (24+70*MAXN)
static EDGE edges[NUMEDGES];

#define init_edge edges

#define quadr_D1(n) (edges -16 + 8*MAXN + 4*(n))
/* The smallest n for which it is called is 3. Then it should start at entry
   edges -4 + 8*MAXN -- right AFTER the edges for quadr_P0(n) which is
   used in the same run. */
#define quadr_P0(n) (edges -4 + 4*MAXN + 4*(n))
/* The smallest n for which it is called is 3. Then it should start at entry
   edges +8 + 4*MAXN -- right AFTER the edges for quadr_P1(n) which is
   used in the same run. */
#define quadr_P1(n) (edges +8 + 4*(n))
/* edges + 24 + 4*(n) - 4*4 since the smallest n for which it is possibly
   called is n=4  and then it should start at edge 24 */

static int degree[MAXN];   /* the degrees of the vertices */
static int class[MAXN]; /* the bipartition class of the vertex */
static int region[MAXN];
static int nregions;
static EDGE *firstedge[MAXN]; /* pointer to arbitrary edge out of vertex i. */
/* This pointer may change during the run, so all one can rely on is that
     at any point it is "some" edge out of i */

static EDGE *numbering[2*MAXE][MAXE];
/* holds numberings produced by canon() or canon_edge() */
static EDGE *saved_numbering[2*MAXE][MAXE];
/* a copy of numbering used by scanordloops() */
static int vindex[MAXN];
/* holds the index of an edge with start i filled in determine_edge_number */
static int adjacency[MAXN][MAXN];
/*The adjacency matrix of the current quadrangulation. The number adjacency[i][j] is the number of edges between i and j*/

static uint16_t edgecode[MAXE + MAXN];
/* precomputed edgecode string that can be output*/
static uint16_t edgecode_mirror[MAXE + MAXN];
/* precomputed edgecode string of the mirrored graph, that can be output*/
static int edgecode_length;
/* The length of edgecode and edgecode_mirror are the same*/


#define PRINTBIG(file,big) fprintf(file,"%llu",(big))


static int markvalue = 30000;
#define RESETMARKS {int mki; if ((markvalue += 2) > 30000) \
       { markvalue = 2; for (mki=0;mki<NUMEDGES;++mki) edges[mki].mark=0;}}
#define MARK(e) ((e)->mark = (markvalue))
#define MARKLO(e) ((e)->mark = (markvalue))
#define MARKHI(e) ((e)->mark = (markvalue+1))
#define UNMARK(e) ((e)->mark = (markvalue-1))
#define ISMARKED(e) ((e)->mark >= markvalue)
#define ISMARKEDLO(e) ((e)->mark == markvalue)
#define ISMARKEDHI(e) ((e)->mark > markvalue)

#define ADD_EDGE(a,b) {adjacency[(a)][(b)]++; adjacency[(b)][(a)]++;}
#define REMOVE_EDGE(a,b) {adjacency[(a)][(b)]--; adjacency[(b)][(a)]--;}


/* and the same for vertices */

static int markvalue_v = 30000;
static int marks__v[MAXN];
#define RESETMARKS_V {int mki; if ((++markvalue_v) > 30000) \
       { markvalue_v = 1; for (mki=0;mki<MAXN;++mki) marks__v[mki]=0;}}
#define UNMARK_V(x) (marks__v[x] = 0)
#define ISMARKED_V(x) (marks__v[x] == markvalue_v)
#define MARK_V(x) (marks__v[x] = markvalue_v)


#define MAX(x,y) ((x)<(y) ? (y) : (x))
#define MIN(x,y) ((x)>(y) ? (y) : (x))

/**************************************************************************/
static void
print_help(){
    fprintf(stderr,"This program generates all non-isomorphic lopsp-operations of a given inflation factor.\n"
                   "and outputs them to stdout in the .lopsp format.\n"
                   "For operations that are not mirror-symmetrical, i.e. they have two chiral forms, only one is given unless option o is used\n"
                   "The options the program takes are:\n"
                   "-a: Operations of smaller inflation factors are also generated. This option is incompatible with -m.\n"
                   "-l: Only lsp-operations are generated.\n"
                   "-o: All operations up to orientation-preserving automorphism are generated.\n"
                   "-s: Some statistics of the operations will be given. Note that this can have a considerable impact on the running time\n"
                   "-u: The operations are not written to stdout. Only the counts are given (through stderr).\n"
                   "-c2 or -c3: Only operations that are 2- resp. 3-connected are generated.\n"
                   "-mi/j where i and j are positive integers such that i < j: This is used to split up large calculations into "
                   "j smaller parts. With mi/j the ith of the j parts will be calculated. For small inflation factors not all parts "
                   "contain quadrangulations. This option is incompatible with -a.\n");
    exit(0);
}

/**************************************************************************/

static void
check_graph_correctness()

/* Method used for tests. Checks if the graph is a valid plane graph */

{
    int edgecount = 0;
    EDGE* edge;
    int deg;
    int adjacencies[MAXN];

    for (int i = 0; i < nv; ++i)
    {
        edge = firstedge[i];
        deg = 0;
        for (int j = 0; j < nv; ++j) {adjacencies[j] = 0;}
        do
        {
            if (class[edge->start] + class[edge->end] != 1){
                fprintf(stderr, "incorrect classes\n");
                exit(1);
            }
            deg++;
            edgecount++;
            adjacencies[edge->end]++;
            if(edge->start != i)
            {
                fprintf(stderr, "inconsistent start\n");
                exit(1);
            }
            if(edge->invers->invers != edge)
            {
                fprintf(stderr, "inconsistent inverses\n");
                exit(1);
            }
            if(edge->invers->end != i)
            {
                fprintf(stderr, "inconsistent end\n");
                exit(1);
            }
            if(edge->next->prev != edge)
            {
                fprintf(stderr, "inconsistent next\n");
                exit(1);
            }
            if(edge->prev->next != edge)
            {
                fprintf(stderr, "inconsistent prev\n");
                exit(1);
            }
            edge = edge->next;
        } while (edge != firstedge[i]);
        if(deg != degree[i])
        {
            fprintf(stderr, "degree not correct\n");
            exit(1);
        }
        for (int j = 0; j < nv; ++j) {
            if(adjacency[i][j] != adjacencies[j] || adjacency[j][i] != adjacencies[j]){
                fprintf(stderr, "adjacency matrix not correct\n");
                exit(1);
            }
        }
    }

    if(edgecount != ne)
    {
        fprintf(stderr, "inconsistent edgecount\n");
        exit(1);
    }
    if((ne/2) + 4 != 2*nv)
    {
        fprintf(stderr, "not a plane quadrangulation\n");
        exit(1);
    }
}

/**************************************************************************/
static void
is_quadrangulation()
/* Method used for tests. Checks if the graph is a valid plane quadrangulation */

{
    check_graph_correctness();
    EDGE* edge;
    int face_size;
    for (int i = 0; i < nv; ++i)
    {
        edge = firstedge[i];
        face_size = 0;
        do
        {
            face_size++;
            edge = edge->invers->next;
        } while (edge != firstedge[i]);
        if(face_size != 4)
        {
            fprintf(stderr, "not a quadrangulation\n");
            exit(1);
        }
    }
}

/*************************CANONICITY*CHECKING*******************************************************************/

/****************************************************************************/

static void
testcanon_first_init(EDGE *givenedge, int representation[], const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->next" direction, an automorphism or even a better representation can
   be found. A better representation will be completely constructed and
   returned in "representation".  It works pretty similar to testcanon except
   for obviously necessary changes, so for extensive comments see testcanon */
{
    register EDGE *run;
    register int vertex;
    EDGE *temp;
    EDGE *startedge[MAXN+1];
    int number[MAXN], i;
    int last_number, actual_number;


    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;
    if (givenedge->start != givenedge->end)
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = givenedge;

    while (last_number < nv)
    {
        for (run = temp->next; run != temp; run = run->next)
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                *representation = colour[vertex]; }
            else *representation = number[vertex];
            representation++;
        }
        *representation = 0;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
    {
        for (run = temp->next; run != temp; run = run->next)
        {
            *representation = number[run->end]; representation++;
        }
        *representation = 0;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }
}

/****************************************************************************/

static void
testcanon_first_init_mirror(EDGE *givenedge, int representation[],
                            const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->prev" direction, an automorphism or even a better representation can
   be found. A better representation will be completely constructed and
   returned in "representation".  It works pretty similar to testcanon except
   for obviously necessary changes, so for extensive comments see testcanon */
{
    register EDGE *run;
    register int vertex;
    EDGE *temp;
    EDGE *startedge[MAXN+1];
    int number[MAXN], i;
    int last_number, actual_number;

    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;
    if (givenedge->start != givenedge->end)
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = givenedge;

    while (last_number < nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                *representation = colour[vertex]; }
            else *representation = number[vertex];
            representation++;
        }
        *representation = 0;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        {
            *representation = number[run->end]; representation++;
        }
        *representation = 0;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }
}

/****************************************************************************/

static int
testcanon_init(EDGE *givenedge, int representation[], const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->next" direction, an automorphism or even a better representation can
   be found. A better representation will be completely constructed and
   returned in "representation".  It works pretty similar to testcanon except
   for obviously necessary changes, so for extensive comments see testcanon */
{
    register EDGE *run;
    register int vertex;
    EDGE *temp;
    EDGE *startedge[MAXN+1];
    int number[MAXN], i;
    int better = 0; /* is the representation already better ? */
    int last_number, actual_number;

    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;
    if (givenedge->start != givenedge->end)
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = givenedge;

    while (last_number < nv)
    {
        for (run = temp->next; run != temp; run = run->next)
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                vertex = colour[vertex]; }
            else vertex=number[vertex];
            if (better) *representation = vertex;
            else {
                if (vertex > (*representation)) return 0;
                else if (vertex < (*representation))
                { better = 1; *representation = vertex; }
            }
            representation++;
        }
        if ((*representation) != 0) { better = 1; *representation = 0; }
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
    {
        for (run = temp->next; run != temp; run = run->next)
        { vertex = number[run->end];
            if (better) *representation = vertex;
            else
            {
                if (vertex > (*representation)) return 0;
                if (vertex < (*representation))
                { better = 1; *representation = vertex; }
            }
            representation++;
        }
        if ((*representation) != 0) { better = 1; *representation = 0; }
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    if (better) return 2;
    return 1;
}

/****************************************************************************/

static int
testcanon_mirror_init(EDGE *givenedge, int representation[], const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->prev" direction, an automorphism or even a better representation can
   be found. A better representation will be completely constructed and
   returned in "representation".  It works pretty similar to testcanon except
   for obviously necessary changes, so for extensive comments see testcanon */
{
    EDGE *temp, *run;
    EDGE *startedge[MAXN+1];
    int number[MAXN], i;
    int better = 0; /* is the representation already better ? */
    int last_number, actual_number, vertex;

    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;
    if (givenedge->start != givenedge->end)
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = givenedge;

    while (last_number < nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                vertex = colour[vertex]; }
            else vertex=number[vertex];
            if (better) *representation = vertex;
            else {
                if (vertex > (*representation)) return 0;
                else if (vertex < (*representation))
                { better = 1; *representation = vertex; }
            }
            representation++;
        }
        if ((*representation) != 0) { better = 1; *representation = 0; }
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        { vertex = number[run->end];
            if (better) *representation = vertex;
            else
            {
                if (vertex > (*representation)) return 0;
                if (vertex < (*representation))
                { better = 1; *representation = vertex; }
            }
            representation++;
        }
        if ((*representation) != 0) { better = 1; *representation = 0; }
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    if (better) return 2;
    return 1;
}

/****************************************************************************/

static void
construct_numb(EDGE *givenedge, EDGE *a_numbering[])

/* Starts at givenedge and writes the edges in the well defined order
   into the list.  Works like testcanon. Look there for comments. */
{
    EDGE *temp, **tail, **limit, *run;
    EDGE *startedge[MAXN+1];
    int last_number, actual_number, vertex;

    RESETMARKS_V

    tail = a_numbering;
    limit = a_numbering+ne-1;

    MARK_V(givenedge->start);
    if (givenedge->start != givenedge->end)
    {
        MARK_V(givenedge->end);
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = *tail = givenedge;

    while (last_number < nv)
    {
        for (run = temp->next; run != temp; run = run->next)
            /* this loop marks all edges around temp->origin. */
        { vertex = run->end;
            if (!ISMARKED_V(vertex))
            { startedge[last_number] = run->invers;
                last_number++; MARK_V(vertex); }
            tail++; *tail = run;
        }
        if (tail != limit)
        { tail++;
            *tail = temp = startedge[actual_number];  actual_number++; }
    }

    while (tail != limit)
        /* Now we know that all numbers have been given */
    {
        for (run = temp->next; run != temp; run = run->next)
            /* this loop marks all edges around temp->origin. */
        { tail++; *tail = run; }
        if (tail != limit)
        {
            /* Next vertex to explore: */
            tail++;
            *tail = temp = startedge[actual_number];  actual_number++; }
    }
}

/****************************************************************************/

static void
construct_numb_mirror(EDGE *givenedge, EDGE *a_numbering[])

/* Starts at givenedge and writes the edges in the well defined order
   into the list.  Works like testcanon. Look there for comments.  */
{
    EDGE *temp, **tail, **limit, *run;
    EDGE *startedge[MAXN+1];
    int last_number, actual_number, vertex;

    RESETMARKS_V

    tail = a_numbering; /* The first entry of the numbering list */
    limit = a_numbering+ne-1;  /* Last valid entry of the numbering list */

    MARK_V(givenedge->start);
    if (givenedge->start != givenedge->end)
    {
        MARK_V(givenedge->end);
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = *tail = givenedge;

    while (last_number < nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
            /* this loop marks all edges around temp->origin. */
        { vertex = run->end;
            if (!ISMARKED_V(vertex))
            { startedge[last_number] = run->invers;
                last_number++; MARK_V(vertex); }
            tail++; *tail = run;
        }
        if (tail != limit)
        {
            tail++;
            *tail = temp = startedge[actual_number];  actual_number++; }
    }

    while (tail != limit)
        /* Now we know that all numbers have been given */
    {
        for (run = temp->prev; run != temp; run = run->prev)
            /* this loop marks all edges around temp->origin. */
        { tail++; *tail = run; }
        if (tail != limit)
        {
            /* Next vertex to explore: */
            tail++;
            *tail = temp = startedge[actual_number];  actual_number++; }
    }
}

/**************************************************************************/

static int
testcanon(EDGE *givenedge, int representation[], const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->next" direction, an automorphism or even a better representation
   can be found. Returns 0 for failure, 1 for an automorphism and 2 for
   a better representation.  This function exits as soon as a better
   representation is found. A function that computes and returns the
   complete better representation can work pretty similar.*/
{
    EDGE *temp, *run;
    EDGE *startedge[MAXN+1]; /* startedge[i] is the starting edge for
                        exploring the vertex with the number i+1 */
    int number[MAXN], i;   /* The new numbers of the vertices, starting
                        at 1 in order to have "0" as a possibility to
                        mark ends of lists and not yet given numbers */
    int last_number, actual_number, vertex;

    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;

    if (givenedge->start != givenedge->end) /* no loop start */
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1 ;

    actual_number = 1;
    temp = givenedge;

/* A representation is a clockwise ordering of all neighbours ended with a 0.
   The neighbours are numbered in the order that they are reached by the BFS
   procedure. In case a vertex is reached for the first time, not the (new)
   number of the vertex is listed, but its colour. When the list around a
   vertex is finished, it is ended with a 0. Since the colours can be
   distinguished from the vertices (requirement for the colour function), the
   adjacency list can be reconstructed: Every time a colour is listed, its
   number would be the smallest number not given yet.
   Since the edges when a vertex is reached for the first time are remembered,
   for these edges we in fact have more than just the vertex information -- for
   these edges we also have the exact information which edge occurs in the
   cyclic order. This makes the routine work also for double edges.

   Since every representation starts with the colour of vertex 2, which is
   the same colour all the time, we do not have to store that.

   In case of a loop as the starting point, the colour of 1 is omitted.
   Nevertheless also in this case it cannot be mixed up with a non-loop
   starting point, since the number of times a colour occurs is larger
   for loop starters than for non-loop starters.

   Every first entry in a new clockwise ordering (the starting point of the
   edge it was numbered from is determined by the entries before (the first
   time it occurs in the list to be exact), so this is not given either.
   The K4 could be denoted  c3 c4 0 4 3 0 2 3 0 3 2 0 if ci is the colour
   of vertex i.  Note that the colour of vertex 1 is -- by definition --
   always the smallest one */

    while (last_number < nv)
    {
        for (run = temp->next; run != temp; run = run->next)
            /* this loop marks all edges around temp->origin. */
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                vertex = colour[vertex]; }
            else vertex = number[vertex];
            if (vertex > (*representation)) return 0;
            if (vertex < (*representation)) return 2;
            representation++;
        }
        /* check whether representation[] is also at the end of a cyclic list */
        if ((*representation) != 0) return 2;
        representation++;
        /* Next vertex to explore: */
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
        /* Now we know that all numbers have been given */
    {
        for (run = temp->next; run != temp; run = run->next)
            /* this loop marks all edges around temp->origin. */
        {
            vertex = number[run->end];
            if (vertex > (*representation)) return 0;
            if (vertex < (*representation)) return 2;
            representation++;
        }
        /* check whether representation[] is also at the end of a cyclic list */
        if ((*representation) != 0) return 2;
        representation++;
        /* Next vertex to explore: */
        temp = startedge[actual_number];  actual_number++;
    }

    return 1;
}

/*****************************************************************************/

static int
testcanon_mirror(EDGE *givenedge, int representation[], const int colour[])

/* Tests whether starting from a given edge and constructing the code in
   "->prev" direction, an automorphism or even a better representation can
   be found. Comments see testcanon -- it is exactly the same except for
   the orientation */
{
    EDGE *temp, *run;
    EDGE *startedge[MAXN+1];
    int number[MAXN], i;
    int last_number, actual_number, vertex;

    for (i = 0; i < nv; i++) number[i] = 0;

    number[givenedge->start] = 1;

    if (givenedge->start != givenedge->end)
    {
        number[givenedge->end] = 2;
        last_number = 2;
        startedge[1] = givenedge->invers;
    }
    else last_number = 1;

    actual_number = 1;
    temp = givenedge;

    while (last_number < nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        { vertex = run->end;
            if (!number[vertex])
            { startedge[last_number] = run->invers;
                last_number++; number[vertex] = last_number;
                vertex = colour[vertex]; }
            else vertex = number[vertex];
            if (vertex > (*representation)) return 0;
            if (vertex < (*representation)) return 2;
            representation++;
        }
        if ((*representation) != 0) return 2;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    while (actual_number <= nv)
    {
        for (run = temp->prev; run != temp; run = run->prev)
        {
            vertex = number[run->end];
            if (vertex > (*representation)) return 0;
            if (vertex < (*representation)) return 2;
            representation++;
        }
        if ((*representation) != 0) return 2;
        representation++;
        temp = startedge[actual_number];  actual_number++;
    }

    return 1;
}


/****************************************************************************/

static int
canon(const int lcolour[], EDGE *can_numberings[][MAXE], int *nbtot, int *nbop)

/* Checks whether the last vertex (number: nv-1) is canonical or not.
   Returns 1 if yes, 0 if not. One of the criterions a canonical vertex
   must fulfill is that its colour is minimal.

   IT IS ASSUMED that the values of the colour function are positive
   and at most INT_MAX-MAXN.

   A possible starting edge for the construction of a representation is
   one with lexicographically minimal colour pair (start,INT_MAX-end).
   In can_numberings[][] the canonical numberings are stored as sequences
   of oriented edges.  For every 0 <= i,j < *nbtot and every
   0 <= k < ne the edges can_numberings[i][k] and can_numberings[j][k] can
   be mapped onto each other by an automorphism. The first
   *nbop numberings are orientation preserving while
   the rest is orientation reversing.

   In case of only 1 automorphism, in can_numberings[0][0] the "canonical"
   edge is given.  It is one edge emanating at the canonical vertex. The
   rest of the numbering is not given.

   In case of nontrivial automorphisms, can[0] starts with a list of edges
   adjacent to nv-1. In case of an orientation preserving numbering the deges
   are listed in ->next direction, otherwise in ->prev direction.

   Works OK if at least one vertex has valence >= 3. Otherwise some numberings
   are computed twice, since changing the orientation (the cyclic order around
   each vertex) doesn't change anything */
{
    int i, last_vertex, test;
    int minstart, maxend; /* (minstart,maxend) will be the chosen colour
                                pair of an edge */
    EDGE *startlist_last[5], *startlist[5*MAXN], *run, *end;
    int list_length_last, list_length;
    int representation[MAXE];
    EDGE *numblist[MAXE], *numblist_mirror[MAXE]; /* lists of edges where
                        starting gives a canonical representation */
    int numbs = 1, numbs_mirror = 0;
    int colour[MAXN];

    for (i=0; i<nv; i++) colour[i]=lcolour[i]+MAXN;
    /* to distinguish colours from vertices */
    last_vertex = nv-1;
    minstart = colour[last_vertex];

/* determine the smallest possible end for the vertex "last_vertex" */

    list_length_last = 1; startlist_last[0] = end = firstedge[last_vertex];
    maxend = colour[end->end];

    for (run = end->next; run != end; run = run->next){
        if (colour[run->end] > maxend){
            startlist_last[0] = run;
            list_length_last = 1; maxend = colour[run->end];
        }
        else if (colour[run->end] == maxend){
            startlist_last[list_length_last] = run; list_length_last++;
        }
    }

/* Now we know the pair that SHOULD be minimal and we can determine a list
   of all edges with this colour pair. If a new pair is found that is even
   smaller, we can return 0 at once */

    list_length = 0;

    for (i = 0; i < last_vertex; i++)
    { if (colour[i] < minstart) return 0;
        if (colour[i] == minstart)
        { run = end = firstedge[i];
            do
            {
                if (colour[run->end] > maxend) return 0;
                if (colour[run->end] == maxend)
                { startlist[list_length] = run; list_length++; }
                run = run->next;
            } while (run != end);
        }
    }

/* OK -- so there is no smaller pair and now we have to determine the
   smallest representation around vertex "last_vertex": */

    testcanon_first_init(startlist_last[0], representation, colour);
    numblist[0] = startlist_last[0];
    test = testcanon_mirror_init(startlist_last[0], representation, colour);
    if (test == 1)
    { numbs_mirror = 1; numblist_mirror[0] = startlist_last[0]; }
    else if (test == 2)
    { numbs_mirror = 1; numbs = 0;
        numblist_mirror[0] = startlist_last[0]; }

    for (i = 1; i < list_length_last; i++)
    { test = testcanon_init(startlist_last[i], representation, colour);
        if (test == 1) { numblist[numbs] = startlist_last[i]; numbs++; }
        else if (test == 2)
        { numbs_mirror = 0; numbs = 1; numblist[0] = startlist_last[i]; }
        test = testcanon_mirror_init(startlist_last[i],
                                     representation, colour);
        if (test == 1)
        { numblist_mirror[numbs_mirror] = startlist_last[i];
            numbs_mirror++; }
        else if (test == 2)
        { numbs_mirror = 1; numbs = 0;
            numblist_mirror[0] = startlist_last[i]; }
    }

/* Now we know the best representation we can obtain starting at last_vertex.
   Now we will check all the others. We can return 0 at once if we find a
   better one */

    for (i = 0; i < list_length; i++)
    { test = testcanon(startlist[i], representation, colour);
        if (test == 1) { numblist[numbs] = startlist[i]; numbs++; }
        else if (test == 2) return 0;
        test = testcanon_mirror(startlist[i], representation, colour);
        if (test == 1)
        { numblist_mirror[numbs_mirror] = startlist[i]; numbs_mirror++; }
        else if (test == 2) return 0;
    }

    *nbop = numbs;
    *nbtot = numbs+numbs_mirror;

    if (*nbtot>1)
    { for (i = 0; i < numbs; i++)
            construct_numb(numblist[i], can_numberings[i]);
        for (i = 0; i < numbs_mirror; i++, numbs++)
            construct_numb_mirror(numblist_mirror[i],can_numberings[numbs]);
    }
    else
    {
        if (numbs) can_numberings[0][0] = numblist[0];
        else can_numberings[0][0] = numblist_mirror[0]; }

    return 1;
}

/****************************************************************************/

static int
canon_edge_oriented(EDGE *edgelist_or[], int num_edges_or, int can_edges_or,
                    EDGE *edgelist_inv[], int num_edges_inv, int can_edges_inv,
                    const int lcolour[], EDGE *can_numberings[][MAXE],
                    int *nbtot, int *nbop)

/*
   IT IS ASSUMED that the values of the colour function are positive
   and at most INT_MAX-MAXN.

   This routine checks all num_edges_or elements of edgelist_or just for one
   orientation and all num_edges_inv elements of the list edgelist_inv just
   for the other. It returns 1 if and only if one of the first can_edges_or
   elements of the first list or first can_edges_inv elements of the second
   give an optimal numbering among all the possibilities provided by the
   lists.

   Edges given are not in minimal form -- but it is guaranteed that all
   colours of the startpoints are the same and all colours of the endpoints
   are the same.

   In case of only the identity automorphism, the entries of can_numberings[][]
   are undefined.

   Otherwise in can_numberings[][] the canonical numberings are stored as
   sequences of oriented edges.  For every 0 <= i,j < *nbtot
   and every 0 <= k < ne the edges can_numberings[i][k] and
   can_numberings[j][k] can be mapped onto each other by an automorphism.
   The first *nbop numberings are orientation
   preserving while the rest are orientation reversing.

   In case of an orientation preserving numbering the edges are listed in
   ->next direction, otherwise in ->prev direction.

   Works OK if at least one vertex has valence >= 3. Otherwise some numberings
   are computed twice, since changing the orientation (the cyclic order around
   each vertex) doesn't change anything */
{
    int i, test;
    int representation[MAXE];
    EDGE *numblist[MAXE], *numblist_mirror[MAXE]; /* lists of edges where
                            starting gives a canonical representation */
    int numbs = 1, numbs_mirror = 0;
    int colour[MAXN];

    for (i=0; i<nv; i++) colour[i]=lcolour[i]+MAXN;
    /* to distinguish colours from vertices */

/* First we have to determine the smallest representation possible with
   edgelist_or */

    if (can_edges_or > 0)
    { testcanon_first_init(edgelist_or[0], representation, colour);
        numblist[0] = edgelist_or[0];
        for (i=1; i<can_edges_or; i++)
        { test = testcanon_init(edgelist_or[i], representation, colour);
            if (test == 1) { numblist[numbs] = edgelist_or[i]; numbs++; }
            else if (test == 2)
            { numbs = 1; numblist[0] = edgelist_or[i]; }
        }
        i=0; /* the next for-loop can start at the beginning */
    }
    else
    { numbs=0; numbs_mirror=1;
        testcanon_first_init_mirror(edgelist_inv[0], representation, colour);
        numblist_mirror[0] = edgelist_inv[0];
        i=1; /* the next for-loop must start at position 1 */
    }

    for (   ; i<can_edges_inv; i++)
    { test = testcanon_mirror_init(edgelist_inv[i], representation, colour);
        if (test == 1)
        { numblist_mirror[numbs_mirror] = edgelist_inv[i]; numbs_mirror++; }
        else if (test == 2)
        { numbs = 0; numbs_mirror=1; numblist_mirror[0] = edgelist_inv[i]; }
    }


    /* now we know the best we can get for a "canonical edge".
       Next we will check all the others. We can return 0 at once if we find a
       better one */

    for (i=can_edges_or ; i < num_edges_or; i++)
    { test = testcanon(edgelist_or[i], representation, colour);
        if (test == 1) { numblist[numbs] = edgelist_or[i]; numbs++; }
        else if (test == 2) return 0;
    }
    for (i=can_edges_inv ; i < num_edges_inv; i++)
    { test = testcanon_mirror(edgelist_inv[i], representation, colour);
        if (test == 1)
        { numblist_mirror[numbs_mirror] = edgelist_inv[i]; numbs_mirror++; }
        else if (test == 2) return 0;
    }

    *nbop = numbs;
    *nbtot = numbs+numbs_mirror;

    if (*nbtot > 1)
    {
        for (i = 0; i < numbs; i++)
            construct_numb(numblist[i], can_numberings[i]);
        for (i = 0; i < numbs_mirror; i++, numbs++)
            construct_numb_mirror(numblist_mirror[i],can_numberings[numbs]);
    }

    return 1;
}

/**************************************************************************/

static void
prune_oriented_lists(EDGE *good_or[], int *ngood_or, int *ngood_ref,
                     EDGE *good_mir[], int *ngood_mir, int *ngood_mir_ref)

/* Try to prune the edge lists (of the form required by
   canon_edge_oriented()) by using the degrees of a couple of
   extra vertices.  The result is returned in the same form.
   As always, if *ngood_ref==*ngood_mir_ref==0 on output
   (which must not be true on input), all else is undefined.
*/
{
    int i,k,kref;
    long code_or[MAXE],code_mir[MAXE],maxcode;
#define PRUNE_OR(e) ((degree[(e)->invers->prev->prev->end] << 10) \
                      + degree[(e)->next->invers->prev->end])
#define PRUNE_MIR(e) ((degree[(e)->invers->next->next->end] << 10) \
                      + degree[(e)->prev->invers->next->end])

    maxcode = 0;
    for (i = 0; i < *ngood_or; ++i)
    {
        code_or[i] = PRUNE_OR(good_or[i]);
        if (code_or[i] > maxcode) maxcode = code_or[i];
    }

    for (i = 0; i < *ngood_mir; ++i)
    {
        code_mir[i] = PRUNE_MIR(good_mir[i]);
        if (code_mir[i] > maxcode) maxcode = code_mir[i];
    }

    k = kref = 0;
    for (i = 0; i < *ngood_or; ++i)
        if (code_or[i] == maxcode)
        {
            if (i < *ngood_ref) ++kref;
            good_or[k++] = good_or[i];
        }
    *ngood_or = k;
    *ngood_ref = kref;

    k = kref = 0;
    for (i = 0; i < *ngood_mir; ++i)
        if (code_mir[i] == maxcode)
        {
            if (i < *ngood_mir_ref) ++kref;
            good_mir[k++] = good_mir[i];
        }
    *ngood_mir = k;
    *ngood_mir_ref = kref;
}

/****************************************************************************/

static void
determine_edge_numbers()

/* Fills in the array vindex such that numbering[0][vindex[i]]->start == i
 * Only call this if nbtot > 1 */

{
    EDGE **nb0;

    int i;

    nb0 = (EDGE**) numbering[0];

    for (i = 0; i < ne; ++i) {
        vindex[nb0[i]->start] = i;
    }
}

///****************************************************************************/
//
//static int
//numoporbits(int nbtot, int nbop)
//
///* return number of orbits of vertices, under the orientation-preserving
//   group (assumed computed). */
//
//{
//    EDGE **nb0,**nb,**nboplim;
//    int vindex[MAXN];
//    int i,count;
//
//    if (nbtot == 1) return nv;
//
//    nb0 = (EDGE**) numbering[0];
//    nboplim = (EDGE**)numbering[nbop==0?nbtot:nbop];
//
//    for (i = 0; i < ne; ++i) vindex[nb0[i]->start] = i;
//
//    RESETMARKS_V
//
//    count = 0;
//    for (i = 0; i < nv; ++i)
//        if (!ISMARKED_V(i))
//        {
//            ++count;
//            for (nb = nb0+vindex[i]; nb < nboplim; nb += MAXE)
//                MARK_V((*nb)->start);
//        }
//
//    return count;
//}

/****************************************************************************/

static void
make_P3(void)

/* Make the path P3 */
{
    EDGE *buffer;

    buffer=edges; /* edge number 0 */
    buffer->start=0;
    buffer->end=1;
    buffer->next=buffer;
    buffer->prev=buffer;
    buffer->invers=edges+1;
    buffer->min=buffer;

    buffer++; /* edge number 1 */
    buffer->start=1;
    buffer->end=0;
    buffer->next=edges+2;
    buffer->prev=edges+2;
    buffer->invers=edges;
    buffer->min=edges;

    buffer++; /* edge number 2 */
    buffer->start=1;
    buffer->end=2;
    buffer->next=edges+1;
    buffer->prev=edges+1;
    buffer->invers=edges+3;
    buffer->min=buffer;

    buffer++; /* edge number 3 */
    buffer->start=2;
    buffer->end=1;
    buffer->next=buffer;
    buffer->prev=buffer;
    buffer->invers=edges+2;
    buffer->min=edges+2;

    firstedge[0]=edges;
    firstedge[1]=edges+1;
    firstedge[2]=edges+3;

    degree[0]=degree[2]=1;
    degree[1]=2;
    class[0]=class[2]=0;
    class[1]=1;

    nv=3;
    ne=4;

    ADD_EDGE(0,1);
    ADD_EDGE(1,2);
}

/**************************************************************************/
static void
initialize_multiquadrangulations(void)
/* initialize edges for the generation of all quadrangulations and
   make the initial path P3 */
{
    EDGE *run,*start;
    int n;

/* Initialize the edges for operation P0. They look like

       a
      / \
    ?b   c    Vertex c is the new point. (a,c) and (d,c) are the
      \ /     new edges with a-->c taken as quadr_P0(n)
       d

It is assumed that for 3<=n<MAXN after quadr_P0(n) there is room for 4 edges.
*/

    for (n=3; n<MAXN; n++)
    { run=start=quadr_P0(n);
        run->end=n; run->min=run;
        run->invers=start+1;

        run=start+1; run->start=n;
        run->prev=run->next=start+3;
        run->min=run->invers=start;

        run=start+2; run->end=n; run->min=run; run->invers=start+3;

        run=start+3; run->start=n; run->min=run->invers=start+2;
        run->prev=run->next=start+1;
    }

/* Then initialize the edges for operation P1. They look like

       a
    ? / \
   --b   c--  Vertex c is the new point. (a,c) and (d,c) are the
    ? \ /     new edges with a-->c taken as quadr_P1(n)
       d

It is assumed that for 5<=n<MAXN after quadr_P1(n) there is room for 4 edges.
*/

    for (n=4; n<MAXN; n++)
    { run=start=quadr_P1(n);
        run->end=n; run->min=run;
        run->invers=start+1;

        run=start+1; run->start=n;
        run->min=run->invers=start; run->prev=start+3;

        run=start+2; run->end=n; run->min=run; run->invers=start+3;

        run=start+3; run->start=n; run->min=run->invers=start+2;
        run->next=start+1;
    }

    /* Finally initialize the edges for operation D1. They look like

           ___
          /   \
       ? a--c  b ?  Vertex c is the new point. (a,b) (upper) and (a,c) are the
          \___/     new edges with a-->c taken as quadr_D1(n)


    It is assumed that for 3<=n<MAXN after quadr_D1(n) there is room for 4 edges.
    */

    for (n=3; n<MAXN; n++)
    { run = start = quadr_D1(n);
        run->end=n; run->min=run;
        run->invers=start+1; run->prev = start + 2;

        run=start+1; run->start=n;
        run->min=run->invers=start; run->prev = run->next = start+1;

        run=start+2; run->min=run; run->invers=start+3; run->next = start;

        run=start+3; run->min=run->invers=start+2;
    }

    make_P3();
}

/*******************************************************************/
static EDGE
*extend_quadr_D1(EDGE *e1)

/* extends a graph in the way given by the type D1 extension for
   quadrangulations.

   The new edge is inserted on the left of e1.

   It returns the new edge ending at the new vertex of degree 1.

   In the picture: e1=a->b below, the edge a->c is returned and
   vertex c is the new point nv.


       ___
      /   \
    ?a--c  b?
      \___/
        e1
*/

{
    register EDGE *start, *dummy;
    int buffer = e1->start;

    start=quadr_D1(nv);

    degree[nv] = 1;
    class[nv] = class[e1->end];

    dummy = e1->prev;
    start->start = buffer;

    start->next = e1; e1->prev = start;
    degree[buffer] += 2;

    //start + 1
    start++;
    firstedge[nv] = start;
    start->end = buffer;

    //start + 2
    start++;
    start->start = buffer; start->end = e1->end;
    start->prev = dummy; dummy->next = start;
    degree[e1->end]++;

    //start + 3
    start++;
    dummy = e1->invers;
    start->end = buffer; start->start = e1->end;
    dummy->next->prev = start; start->next = dummy->next;
    start->prev = dummy; dummy->next = start;

    ADD_EDGE(nv,buffer);
    ADD_EDGE(buffer,e1->end);

    nv++; ne+=4;

    return (start - 3);
}
/**************************************************************************/

static void
reduce_quadr_D1(EDGE *e)

/* reduces a graph previously extended by the D1 extension for
   quadrangulations. The edge e must be the reference edge
   returned by the expansion routine.
*/

{
    register EDGE *dummy, *dummy2;
    int buffer;

    nv--; ne-=4;

    dummy = e->next; //dummy is the lower edge a->b
    dummy2 = e->prev->prev;

    dummy2->next = dummy; dummy->prev = dummy2;
    buffer = dummy->start;
    firstedge[buffer] = dummy; degree[buffer] -= 2;

    dummy = dummy->invers;
    dummy2 = dummy->next->next;

    dummy->next = dummy2; dummy2->prev = dummy;
    buffer = dummy->start;
    firstedge[buffer] = dummy; degree[buffer]--;

    REMOVE_EDGE(e->start, e->end);
    REMOVE_EDGE(e->start, buffer);
}
/*******************************************************************/

static EDGE
*extend_quadr_P0(EDGE *e1)

/* extends a graph in the way given by the type P0 extension for
   quadrangulations.

   The new path is inserted on the left of e1.

   It returns the directed edge starting at the new vertex with the
   new face on its left.

   In the picture: e1=b->a, the edge c->a is returned and
   vertex c is the new point nv.

       ?
       a
      / \
    ?b   c?
      \ /
       d
       ?
*/

{
    register EDGE *start, *dummy, *dummy2;
    int buffer;

    start=quadr_P0(nv);

    degree[nv]=2;
    class[nv] = class[e1->start];

    buffer=e1->end;
    dummy=e1->invers; dummy2=dummy->prev;
    start->start=buffer;
    start->next=dummy; dummy->prev=start;
    start->prev=dummy2; dummy2->next=start;
    degree[buffer]++;
    ADD_EDGE(nv,buffer);

    start++;
    firstedge[nv]=start;
    start->end=buffer;

    e1=e1->next;
    start++;
    buffer=e1->end;
    dummy=e1->invers; dummy2=dummy->next;
    start->start=buffer;
    start->next=dummy2; dummy2->prev=start;
    start->prev=dummy; dummy->next=start;
    degree[buffer]++;
    ADD_EDGE(nv,buffer);

    (start+1)->end=buffer;

    nv++; ne+=4;

    return (start-1);
}

/**************************************************************************/

static void
reduce_quadr_P0(EDGE *e)

/* reduces a graph previously extended by the P0 extension for
   quadrangulations. The edge e must be the reference edge
   returned by the expansion routine.
*/

{
    register EDGE *dummy, *dummy2;
    int buffer;

    nv--; ne-=4;

    dummy=e->invers;
    dummy2=dummy->next;
    dummy=dummy->prev;

    dummy->next=dummy2; dummy2->prev=dummy;
    buffer=dummy->start;
    firstedge[buffer]=dummy; degree[buffer]--;
    REMOVE_EDGE(buffer, nv)

    dummy=e->prev->invers;
    dummy2=dummy->next;
    dummy=dummy->prev;

    dummy->next=dummy2; dummy2->prev=dummy;
    buffer=dummy->start;
    firstedge[buffer]=dummy; degree[buffer]--;
    REMOVE_EDGE(buffer, nv)
}

/*******************************************************************/

static EDGE
*extend_quadr_P1(EDGE *e)

/* extends a graph in the way given by the type P1 extension for
   3-connected quadrangulations.

   It returns the directed edge characterizing this operation.
*/

{
    register EDGE *e1, *e1i, *e2, *e3, *e3i, *e4, *start, *dummy;
    int end1, end2, center;

    REMOVE_EDGE(e->start, e->end);
    ADD_EDGE(nv,e->end);
    class[nv] = class[e->start];

    center=e->start;
    e1=e->prev;
    e1i=e1->invers;
    e2=e1i->prev;
    end1=e2->start;

    e3=e->next;
    e3i=e3->invers;
    e4=e3i->next;
    end2=e4->start;

    e1->next=e3; e3->prev=e1; degree[center]--; firstedge[center]=e1;

    start=quadr_P1(nv);
    dummy=start+1;
    start->start=dummy->end=end1;
    degree[end1]++;
    start->next=e1i; e1i->prev=start;
    start->prev=e2; e2->next=start;
    ADD_EDGE(nv,end1);

    dummy->next=e; e->prev=dummy;

    start+=2; dummy=start+1;
    start->start=dummy->end=end2;
    degree[end2]++;
    e3i->next=start; start->prev=e3i;
    e4->prev=start; start->next=e4;

    e->next=dummy; dummy->prev=e;

    degree[nv]=3; e->start=e->invers->end=nv; firstedge[nv]=e;
    ADD_EDGE(nv,end2);
    nv++; ne+=4;

    return (e);
}

/**************************************************************************/

static void
reduce_quadr_P1(EDGE *e)

/* reduces a graph previously extended by the type P1 extension for
   3-connected quadrangulations. The edge e must be the reference edge
*/

{
    register EDGE *e1, *e1i, *e2, *e3, *e3i, *e4;
    int end1, end2, center;

    nv--; ne-=4;

    REMOVE_EDGE(e->start, e->end);

    e2=e->prev->invers->prev;
    end1=e2->start;
    e1i=e2->next->next;
    e1=e1i->invers;
    center=e1->start;
    e3=e1->next;
    e3i=e3->invers;
    e4=e3i->next->next;
    end2=e4->start;
    ADD_EDGE(center,e->end);
    REMOVE_EDGE(nv, end1);
    REMOVE_EDGE(nv,end2);

    e2->next=e1i; e1i->prev=e2; firstedge[end1]=e2; degree[end1]--;
    e1->next=e; e->prev=e1; e->next=e3; e3->prev=e;
    e->start=e->invers->end=center; degree[center]++;
    e3i->next=e4; e4->prev=e3i; firstedge[end2]=e4; degree[end2]--;
}

/****************************************************************************/
static void
find_extensions_multiquad(int nbtot, int nbop, EDGE *extD1[], int *nextD1, EDGE *extP0[], int *nextP0,
                          EDGE *extP1[], int *nextP1)

/* Determine the inequivalent places to make extensions, for the
   multiquadrangulations.  The results are put in the arrays
   extP1[0..*nextP1-1], etc..
   nbtot and nbop are the usual group parameters.
*/
{
    EDGE *e, *elast;
    EDGE **nb, **nb0, **nblim, **nboplim;
    int i, j, k, l, x, y;
    int deg1, deg2;
#define VCOLP0(i, j) (degree[i]<degree[j] ? \
 (degree[j]<<10)+degree[i] : (degree[i]<<10)+degree[j])

    if (degree[nv - 1] == 2 && firstedge[nv - 1]->end != firstedge[nv - 1]->next->end) {
        x = firstedge[nv - 1]->end;
        y = firstedge[nv - 1]->next->end;
    } else
        x = -1;

    deg1 = 0;
    deg2 = 0;
    int deg1_list[2];

    int max_neighbourdeg = 0;
    int second_neighbourdeg = 0;
    int v1type1 = !aswitch && (factor % 2 == 1);

    for (i = nv - 1; i >= 0 && deg1 < 4; i--)
        //Once four vertices of degree 1 have been found we have enough information.
        // As the canonical D1 always has the neighbour with the largest degree
        // and we traverse the vertices from last added to first added, max_neighbourdegree
        // is also correct
    {
        if (degree[i] == 2 && firstedge[i]->end != firstedge[i]->next->end) {
            ++deg2;
        }
        if (degree[i] == 1) {
            k = degree[firstedge[i]->end];
            if (deg1 < 2){
                if (deg1 == 0){
                    max_neighbourdeg = k;
                } else {
                    second_neighbourdeg = k;
                }
                deg1_list[deg1] = i;
            }
            ++deg1;
        }
    }

    //if there are 3 or more vertices of degree 1, that number never decreases, so these quadrangulations must not be generated
    if (min_connectivity >= 2){
        if (deg1 > 3 || (deg1 > 2 && min_connectivity == 3)){
            *nextD1 = *nextP0 = *nextP1 = 0;
            return;
        }
    }
    //if v1type1, then we always do a D1 extension in the last step. During the rest of the generation, deg1 can be kept lower

    if (nbtot == 1) {
        /* D1 extension for trivial group */

        k = 0;
        for (l = 0; l < nv; ++l) {
            e = elast = firstedge[l];
            do {
                if (degree[l] == 1){
                    //maxdeg may not be correct after extension
                    //This is also the case if degree[e->end] == 1, but then the
                    //else statement will give the correct answer
                    if (second_neighbourdeg <= 3 && degree[e->end] == max_neighbourdeg){
                        //The degree of the new neighbour will be 3. If there is a larger one e is not canonical
                        extD1[k++] = e;
                    }
                } else {//maxdeg will remain the same or larger after extension
                    if (degree[l] + 2 >= max_neighbourdeg) {
                        extD1[k++] = e;
                    }
                    //we do not look at the degree of the opposite vertex here
                    //as it is difficult to determine which degrees are allowed
                    //(the degree of the neighbour vertex can increase)
                }
                e = e->next;
            } while (e != elast);
        }

        *nextD1 = k;

        /* P0 extension for trivial group */
        //For every graph that is not P3, P0 can remove at most one vertex of degree 1
        if ((deg1 > 1 && nv != 3) || (v1type1 && nv == maxnv - 1))
            *nextP0 = 0;
        else {
            k = 0;
            if (deg1 == 1){
                //there are two possible edges that can be chosen. However, they extend to the same
                //graph, with a different marked edge. One of them always has a lower degree of ->end
                //so we pick the other one
                e = firstedge[deg1_list[0]]->invers->prev;
                i = e->end;
                j = e->next->end;
                //if the last operation applied was not P0 or if the edge was made in the last operations
                if (x < 0 || i == nv - 1 || j == nv - 1)
                    extP0[k++] = e;
                else { //if the last operation applies was P1 then we can check if the new one can be the canonical one or if the previous operation will always be canonical
                    ++degree[i];
                    ++degree[j];
                    if (VCOLP0(i, j) >= VCOLP0(x, y)) extP0[k++] = e;
                    --degree[i];
                    --degree[j];
                }
            } else {
                RESETMARKS
                for (l = 0; l < nv; ++l) {
                    e = elast = firstedge[l];
                    do {
                        if (!ISMARKEDLO(e)) {
                            i = e->end;
                            j = e->next->end;
                            if (i != j) {//not parallel edges
                                //if the last operation applied was not P0 or if the edge was made in the last operations
                                if (x < 0 || i == nv - 1 || j == nv - 1)
                                    extP0[k++] = e;
                                else {
                                    //if the last operation applied was P1 then we must check if the new one can be the
                                    // canonical one or if the previous operation will always be canonical
                                    ++degree[i];
                                    ++degree[j];
                                    if (VCOLP0(i, j) >= VCOLP0(x, y)) extP0[k++] = e;
                                    --degree[i];
                                    --degree[j];
                                }
                            }
                            MARKLO(e->next->invers->next->invers);
                        }
                        e = e->next;
                    } while (e != elast);
                }
            }
            *nextP0 = k;
        }
        /* P1 extension for trivial group */
        //P1 can remove at most two vertices of degree 1 and of degree 2 with different neighbours
        // so if there are three the resulting graph can be constructed with D1 or P0
        if ((deg1 + deg2 > 2) || (v1type1 && nv == maxnv - 1))
            *nextP1 = 0;
        else
        {
            k = 0;
            if (deg1 == 2){//two vertices of degree 1 must be removed
                e = firstedge[deg1_list[0]]->invers;
                elast = firstedge[deg1_list[1]]->invers;
                if(e->start == elast->start && degree[e->start] >= 4){
                    if(e->next->next == elast){
                        extP1[k++] = e->next;
                    }
                    if(elast->next->next == e){
                        extP1[k++] = elast->next;
                    }
                } //else not possible to remove both
            } else if (deg1 == 1){//one vertex of degree 1 must be removed
                e = firstedge[deg1_list[0]]->invers;
                if(degree[e->start] >= 4){
                    extP1[k++] = e->next;
                    extP1[k++] = e->prev;
                }
            } else {//no degree 1
                for (l = 0; l < nv; ++l)
                    if (degree[l] >= 4) {
                        e = elast = firstedge[l];
                        do {
                            i = e->next->end;
                            j = e->prev->end;
                            //no parallel edges in face
                            if (i != j) {
                                extP1[k++] = e;
                            }
                            e = e->next;
                        } while (e != elast);
                    }
            }
            *nextP1 = k;
        }
    }
    else
    {
        nboplim = (EDGE**)numbering[nbop==0?nbtot:nbop];
        nblim = (EDGE**)numbering[nbtot];
        nb0 = (EDGE**)numbering[0];

        for (i = 0; i < ne; ++i) nb0[i]->index = i;

        /* D1 extensions for non-trivial group */

        k = 0;
        RESETMARKS

        for (l = 0; l < ne; ++l)
            if (!ISMARKED(nb0[l]))
            {
                e = nb0[l];
                if (degree[e->start] == 1){
                    if (second_neighbourdeg <= 3 && degree[e->end] == max_neighbourdeg){
                        extD1[k++] = e;
                    }
                } else {
                    if (degree[e->start] + 2 >= max_neighbourdeg) {
                        extD1[k++] = e;
                    }
                }
                //no need to mark the edge itself
                //mark all edges in the orientation-preserving orbit
                for (nb = nb0+l+MAXE; nb < nboplim; nb += MAXE)
                    MARKLO(*nb);
                //mark all edges in the orientation-reversing orbit
                for ( ; nb < nblim; nb += MAXE)
                    MARKLO((*nb));
            }

        *nextD1 = k;

        /* P0 extensions for non-trivial group */
        if ((deg1 > 1 && nv != 3) || (v1type1 && nv == maxnv - 1))
            *nextP0 = 0;
            //Testing for deg1 == 1 is almost useless here. In that case there are almost certainly no automorphisms
        else {
            k = 0;
            RESETMARKS

            for (l = 0; l < ne; ++l)
                if (!ISMARKED(nb0[l]))
                {
                    e = nb0[l];
                    i = e->end;
                    j = e->next->end;
                    if(i != j) {//not parallel edges
                        if (x < 0 || i == nv - 1 || j == nv - 1)
                            extP0[k++] = e;
                        else {
                            ++degree[i];
                            ++degree[j];
                            if (VCOLP0(i, j) >= VCOLP0(x, y)) extP0[k++] = e;
                            --degree[i];
                            --degree[j];
                        }
                    }
                    for (nb = nb0+l+MAXE; nb < nboplim; nb += MAXE)
                        MARKLO(*nb);
                    for ( ; nb < nblim; nb += MAXE)
                        MARKLO((*nb)->invers->next->invers);
                    for (nb = nb0+(nb0[l]->next->invers->next->invers->index);
                         nb < nboplim; nb += MAXE)
                        MARKLO(*nb);
                    for ( ; nb < nblim; nb += MAXE)
                        MARKLO((*nb)->invers->next->invers);
                }

            *nextP0 = k;
        }
        /* P1 extensions for non-trivial group */

        if ((deg1 + deg2 > 2) || (v1type1 && nv == maxnv - 1))
            *nextP1 = 0;
        else
        {
            k = 0;
            if (deg1 == 2){//two vertices of degree 1 must be removed
                e = firstedge[deg1_list[0]]->invers;
                elast = firstedge[deg1_list[1]]->invers;
                if(e->start == elast->start && degree[e->start] >= 4){
                    if(e->next->next == elast){
                        extP1[k++] = e->next;
                    }
                    if(elast->next->next == e){
                        extP1[k++] = elast->next;
                    }
                } //else not possible to remove both
            } else if (deg1 == 1){//one vertex of degree 1 must be removed
                e = firstedge[deg1_list[0]]->invers;
                if(degree[e->start] >= 4){
                    extP1[k++] = e->next;
                    extP1[k++] = e->prev;
                }
            } else {//no degree 1
                RESETMARKS
                for (l = 0; l < nv; ++l) {
                    if (degree[l] >= 4) {
                        e = elast = firstedge[l];
                        do {
                            if (!ISMARKEDLO(e)) {
                                i = e->next->end;
                                j = e->prev->end;
                                if (i != j) {
                                    extP1[k++] = e;
                                }
                                for (nb = nb0 + e->index + MAXE; nb < nblim; nb += MAXE)
                                    MARKLO(*nb);
                            }
                            e = e->next;
                        } while (e != elast);
                    }
                }
            }
            *nextP1 = k;
        }
    }
}
/****************************************************************************/
static int
has_quadr_D1(void)

/* Test whether there is a legal D1 reduction */

{
    int i;

    for (i = nv - 1; --i >= 0; )
        if (degree[i] == 1) return TRUE;

    return FALSE;
}
/****************************************************************************/
/*in multiquadrangulations not all vertices of degree 2 can be reduced. Only those
    in simple faces. */
static int
has_quadr_P0_multi(void)

/* Test whether there is a legal P0 reduction */

{
    int i;

    for (i = nv; --i >= 0; )
        if (degree[i] == 2 && firstedge[i]->end != firstedge[i]->next->end) return TRUE;

    return FALSE;
}

/****************************************************************************/

static void
quadr_D1_legal_all(EDGE *ref, EDGE *good_or[], int *ngood_or, int *ngood_ref,
                   EDGE *good_mir[], int *ngood_mir, int *ngood_mir_ref)

/* The D1-operation with reference edge *ref has just been performed.
   Make a list in good_or[0..*ngood_or-1] of the reference edges of
   legal D1-reductions (oriented editions) that might be canonical,
   with the first *ngood_ref of those being ref.
   Make a list in good_mir[0..*ngood_mir-1] of the
   reference edges of legal four-reductions (mirror-image editions)
   that might be canonical, with the first *ngood_mir_ref of those being
   ref->next.
   *ngood_ref and *ngood_mir_ref might each be 0-1.  If they are
   both 0, nothing else need be correct.
   All the edges in good_or[] and good_mir[] must start with the same
   vertex degree and end with the same vertex degree (actually, colour
   as passed to canon_edge_oriented).

    Of course an edge can only be canonical if its end has degree 1
    We define that the start of a canonical edge has a degree that is as high as possible.
    If there are multiple edges with the same degree start, the other vertex in the 2-cycle
    must have degree as high as possible
*/

{
    EDGE *e;
    int maxdeg = degree[ref->start];
    int max_oppositedeg = degree[ref->next->end];
    int i,nor,nmir;

    nor = nmir = 0; //temporary count. If there turns out to be a larger degree, final value is 0

    good_or[nor++] = ref;
    good_mir[nmir++] = ref;

    *ngood_ref = 1;
    *ngood_mir_ref = 1;

    for (i = nv-1; --i >= 0; )
    {
        e = firstedge[i]->invers;
        if (degree[i] == 1) //candidate ends in a vertex of degree 1
        {
            if(degree[e->start] > maxdeg) {
                //better candidate, ref is not canonical
                *ngood_ref = 0;
                *ngood_mir_ref = 0;
                return;
            }
            if (degree[e->start] == maxdeg) {
                if (degree[e->next->end] > max_oppositedeg){ //better candidate, ref is not canonical
                    *ngood_ref = 0;
                    *ngood_mir_ref = 0;
                    return;
                }
                if (degree[e->next->end] == max_oppositedeg){
                    good_or[nor++] = e;
                    good_mir[nmir++] = e;
                }
            }
        }
    }
    if (nor > *ngood_ref || nmir > *ngood_mir_ref)
        prune_oriented_lists(good_or,&nor,ngood_ref,
                             good_mir,&nmir,ngood_mir_ref);

    *ngood_or = nor;
    *ngood_mir = nmir;
}

/****************************************************************************/

static void
multiquadr_P0_legal_all(EDGE *ref, EDGE *good_or[], int *ngood_or, int *ngood_ref,
                        EDGE *good_mir[], int *ngood_mir, int *ngood_mir_ref)

/* The P0-operation with reference edge *ref has just been performed.
   Make a list in good_or[0..*ngood_or-1] of the reference edges of
   legal P0-reductions (oriented editions) that might be canonical,
   with the first *ngood_ref of those being ref.
   Make a list in good_mir[0..*ngood_mir-1] of the
   reference edges of legal four-reductions (mirror-image editions)
   that might be canonical, with the first *ngood_mir_ref of those being
   ref->next.
   *ngood_ref and *ngood_mir_ref might each be 0-1.  If they are
   both 0, nothing else need be correct.
   All the edges in good_or[] and good_mir[] must start with the same
   vertex degree and end with the same vertex degree (actually, colour
   as passed to canon_edge_oriented).
   P0-reductions have a priority (colour) based on the degrees of the
   end vertex and two side vertices.  It cannot be changed without
   changing find_extensions_quad too.

   version for general quadrangulations that may not be simple
   if there is a vertex of degree 2 its neighbours must be different for it to be a candidate
*/

{
    EDGE *e;
    int col,col2,maxcol,maxcol2;
    int i,nor,nmir;
    int d1,d2,d3,d4;

#define VCOLPD(di,dj) ((di)<(dj) ? ((dj)<<10)+(di) : ((di)<<10)+(dj))

    nor = nmir = 0;

    d1 = degree[ref->end];
    d2 = degree[ref->next->end];
    d3 = degree[ref->invers->next->end];
    d4 = degree[ref->invers->prev->end];

    maxcol = VCOLPD(d1,d2);
    maxcol2 = VCOLPD(d3,d4);

    if (d1 >= d2) //this edge could be candidate
    {
        if (d3 >= d4) good_or[nor++] = ref;
        if (d4 >= d3) good_mir[nmir++] = ref;
    }
    if (d2 >= d1) // the other edge at the degree 2 vertex is a candidate
    {
        if (d4 >= d3) good_or[nor++] = ref->next;
        if (d3 >= d4) good_mir[nmir++] = ref->next;
    }

    *ngood_ref = nor;
    *ngood_mir_ref = nmir;

    for (i = nv-1; --i >= 0; )
        if (degree[i] == 2 && firstedge[i]->end != firstedge[i]->next->end)
        {
            e = firstedge[i];
            d1 = degree[e->end];
            d2 = degree[e->next->end];
            col = VCOLPD(d1,d2);
            if (col > maxcol)
            {
                *ngood_ref = *ngood_mir_ref = 0;
                return;
            }
            else if (col == maxcol)
            {
                d3 = degree[e->invers->next->end];
                d4 = degree[e->invers->prev->end];
                col2 = VCOLPD(d3,d4);

                if (col2 > maxcol2)
                {
                    *ngood_ref = *ngood_mir_ref = 0;
                    return;
                }
                else if (col2 == maxcol2)
                {
                    if (d1 >= d2)
                    {
                        if (d3 >= d4) good_or[nor++] = e;
                        if (d4 >= d3) good_mir[nmir++] = e;
                    }
                    if (d2 >= d1)
                    {
                        if (d4 >= d3) good_or[nor++] = e->next;
                        if (d3 >= d4) good_mir[nmir++] = e->next;
                    }
                }
            }
        }

    if (nor > *ngood_ref || nmir > *ngood_mir_ref)
        prune_oriented_lists(good_or,&nor,ngood_ref,
                             good_mir,&nmir,ngood_mir_ref);

    *ngood_or = nor;
    *ngood_mir = nmir;
}

/****************************************************************************/

static void
multiquadr_P1_legal_all(EDGE *ref, EDGE *good_or[], int *ngood_or, int *ngood_ref,
                        EDGE *good_mir[], int *ngood_mir, int *ngood_mir_ref)

/* The P1-operation with reference edge *ref has just been performed.
   Make a list in good_or[0..*ngood_or-1] of the reference edges of
   legal P1-reductions (oriented editions) that might be canonical,
   with the first *ngood_ref of those being ref.
   Make a list in good_mir[0..*ngood_mir-1] of the
   reference edges of legal four-reductions (mirror-image editions)
   that might be canonical, with the first *ngood_mir_ref of those being
   ref->next.
   *ngood_ref and *ngood_mir_ref might each be 0-1.  If they are
   both 0, nothing else need be correct.
   All the edges in good_or[] and good_mir[] must start with the same
   vertex degree and end with the same vertex degree (actually, colour
   as passed to canon_edge_oriented).
   P1-reductions have a priority (colour) based on the degrees of the
   end vertex and two side vertices.  It cannot be changed without
   changing find_extensions_quad too.

   version for general quadrangulations
   if there is a vertex of degree 3 two of its neighbours must be different for the other edge to be a candidate
*/

{
    EDGE *e,*ee,*elast;
    int maxdeg;
    int col,maxcol;
    int i,nor,nmir;

#define QORCOL(e) (((long)degree[(e)->prev->end]<<10)+(long)degree[(e)->next->end])
#define QMIRCOL(e) (((long)degree[(e)->next->end]<<10)+(long)degree[(e)->prev->end])
//checks whether e can be the marked edge of a legal P1 reduction. No parallel edges
#define LEGAL_P1(e) (((e)->next->end != (e)->prev->end) && ((e)->start != (e)->next->invers->prev->end))

    RESETMARKS
    nor = nmir = 0;

    maxdeg = degree[ref->end];
    maxcol = QORCOL(ref);
    col = QMIRCOL(ref);
    if (col < maxcol)
        good_or[nor++] = ref;
    else if (col == maxcol)
    {
        good_or[nor++] = ref;
        good_mir[nmir++] = ref;
    }
    else
    {
        maxcol = col;
        good_mir[nmir++] = ref;
    }

    *ngood_ref = nor;
    *ngood_mir_ref = nmir;

    MARKLO(ref->invers);

    if (nor > *ngood_ref || nmir > *ngood_mir_ref)
        prune_oriented_lists(good_or,&nor,ngood_ref,
                             good_mir,&nmir,ngood_mir_ref);

    if (*ngood_ref == 0 && *ngood_mir_ref == 0) return;

    for (i = 0; i < nv; ++i)
        if (degree[i] > maxdeg)
        {
            e = elast = firstedge[i];
            do
            {
                ee = e->invers;
                if (!ISMARKEDLO(e) && degree[e->end] == 3 && LEGAL_P1(ee))
                {
                    *ngood_ref = *ngood_mir_ref = 0;
                    return;
                }
                e = e->next;
            } while (e != elast);
        }
        else if (degree[i] == maxdeg)
        {
            e = elast = firstedge[i];
            do
            {
                ee = e->invers;
                if (!ISMARKEDLO(e) && degree[e->end] == 3 && LEGAL_P1(ee))
                {
                    col = QORCOL(e->invers);
                    if (col > maxcol)
                    {
                        *ngood_ref = *ngood_mir_ref = 0;
                        return;
                    }
                    else if (col == maxcol)
                        good_or[nor++] = e->invers;

                    col = QMIRCOL(e->invers);
                    if (col > maxcol)
                    {
                        *ngood_ref = *ngood_mir_ref = 0;
                        return;
                    }
                    else if (col == maxcol)
                        good_mir[nmir++] = e->invers;
                }
                e = e->next;
            } while (e != elast);
        }

    if (nor > *ngood_ref || nmir > *ngood_mir_ref)
        prune_oriented_lists(good_or,&nor,ngood_ref,
                             good_mir,&nmir,ngood_mir_ref);

    *ngood_or = nor;
    *ngood_mir = nmir;
}

/**************************************************************************/
static EDGE*
next_marked_edge(EDGE* edge){
    //assumes that there is a marked edge with start edge->start
    EDGE* result = edge->next;
    int count = 0;
    while (!ISMARKED(result) && count <= degree[edge->start]) {result = result->next;count++;}
    if (!ISMARKED(result)){
        fprintf(stderr, "dju\n");
    }
    return result;
}

/**************************************************************************/
static int
marked_facesize(EDGE* edge){
    int length = 0;
    EDGE* faceedge = edge;
    do {
        length++;
        faceedge = next_marked_edge(faceedge->invers);
    } while (faceedge != edge);
    return length;
}

/**************************************************************************/
static void
unindex_markedface(EDGE* edge){
    EDGE* faceedge = edge;
    do {
        faceedge->index = 0;
        faceedge = next_marked_edge(faceedge->invers);
    } while (faceedge != edge);
}

/**************************************************************************/
static void
mark_edges_to_marked_vertices(int vertex){
    //marks all edges (and their inverses) starting at vertex that end in a marked vertex
    EDGE* edge = firstedge[vertex];
    do {
        if(ISMARKED_V(edge->end)){
            MARK(edge);
            MARK(edge->invers);
        }
        edge = edge->next;
    } while (edge != firstedge[vertex]);
}


/**************************************************************************/
static int
find_vi_recursive(EDGE* edge, uint16_t v0, uint16_t v1, uint16_t v2) {
    /*The vertices of a face have been marked. edge is an edge on the inside of the facial cycle. Check if one of the
     * vi is inside. This method must be called multiple times if not all vertices can be reached from edge*/
    if (ISMARKED_V(edge->end)){
        return FALSE;
    }
    if (edge->end == v0 || edge->end == v1 || edge->end == v2){
        return TRUE;
    }
    MARK_V(edge->end);
    EDGE* neighbour = edge->invers->next;
    while (neighbour != edge->invers){
        if (find_vi_recursive(neighbour, v0, v1, v2)){
            return TRUE;
        }
        neighbour = neighbour->next;
    }
    return FALSE;
}

/**************************************************************************/
static int
vi_inside_or_empty(EDGE* edge, uint16_t v0, uint16_t v1, uint16_t v2) {
    /*A subgraph is marked. Determine if v0, v1 or v2 are in the face of that subgraph containing edge edge */
    EDGE* faceedge = edge;
    EDGE* inneredge;
    int number_recursive_called = 0;
    int found = FALSE;

    do {
        inneredge = faceedge->prev;

        while (!ISMARKED(inneredge) ){
            if (!ISMARKED_V(inneredge->end)) {//this must be checked as it is possible that there is an edge inside but no vertices. In that case it is empty
                if (found){ //to make sure we do not look at this face again
                    MARK_V(inneredge->end);
                } else if (find_vi_recursive(inneredge, v0, v1, v2)) {
                    found = TRUE;
                }
                number_recursive_called++;
            }
            inneredge = inneredge->prev;
        }
        faceedge = next_marked_edge(faceedge->invers);
    } while (faceedge != edge);
    return found || (number_recursive_called == 0);
}

/**************************************************************************/
static int
assign_region_rec(EDGE* from, int number, int* visited) {
    //a recursive method for marking all vertices in a region.
    //returns number if the region was succesfully assigned this number.
    //returns the number of the region already in place if not, and 0 if this recursion branch found no vertices in the region
    //visited[v] is TRUE if v was already visited by this method
    int result;
    int found = FALSE;
    int vertex = from->end;
    EDGE* edge;
    if (visited[vertex]){
        return -1;
    }
    visited[vertex] = TRUE;
    if (region[vertex] == -1){
        found = TRUE;
        region[vertex] = number;
        edge = firstedge[vertex];
        do {
            if ((result = assign_region_rec(edge, number, visited)) != number){
                if (result != -1){
                    return result;
                }
            }
            edge = edge->next;
        } while (edge != firstedge[vertex]);
    } else if (ISMARKED_V(vertex)){
        // a vertex in the subgraph. Continue but only for edges between marked edges
        edge = from->invers->prev;
        while (!ISMARKED(edge)){edge = edge->prev;}
        edge = edge->next;
        while (!ISMARKED(edge)){
            if ((result = assign_region_rec(edge, number, visited)) != number){
                if (result != -1){
                    return result;
                }
            } else {
                found = TRUE;
            }
            edge = edge->next;
        }
    } else {
        //the vertex should be unmarked but isn't. Therefore this region has been seen before
        return region[vertex];
    }
    return found ? number : -1;
}

/**************************************************************************/
static int
assign_region(EDGE* edge, int number) {
    //returns number if the region was successfully assigned this number.
    //returns the number of the region already in place if not, and 0 if there are no vertices in the region
    EDGE* inneredge, *stopedge;
    int visited[nv];
    for (int i = 0; i < nv; ++i) {visited[i] = 0;}
    int result;
    EDGE* faceedge = edge;
    stopedge = edge->prev;
    int found = FALSE;
    while (!ISMARKED(stopedge->invers)){stopedge = stopedge->prev;}
    do {
        inneredge = faceedge->prev;
        while (inneredge != stopedge){
            if ((result = assign_region_rec(inneredge, number, visited)) != number){
                if (result != -1){
                    return result;
                }
            } else {
                found = TRUE;
            }
            inneredge = inneredge->prev;
        }

        stopedge = faceedge->invers;
        faceedge = next_marked_edge(faceedge->invers);
    } while (faceedge != edge);
    return found ? number : 0;
}

/**************************************************************************/
//Does some easy checks to see if the lopsp-operation is certainly not c3 or maybe c3
static int
maybe_c3(uint16_t v0, uint16_t v1, uint16_t v2, int v1type1) {
    if (nv == 3 && v1type1){//P3 is a special case when it comes to cycles, so it is treated separately
        return TRUE;
    }
    if (adjacency[v0][v2]){//there can be no edges between v0 and v2
        return FALSE;
    }
    if (nv != 3 && adjacency[v0][v1] && adjacency[v1][v2]){
        return FALSE;
    }

    for (int v = 0; v < nv; ++v) {
        if (degree[v] == 2) {
            if(v != v0 && v != v1 && v != v2){
                return FALSE;
            }
        }
        for (int w = v + 1; w < nv; ++w) {
            switch (adjacency[v][w]) {
                case 0:
                case 1:
                    break;
                case 2:
                    if (v == v0 || v == v2  || w == v0 || w == v2){ //one of the faces cannot have v0 or v2 -> v1 must be alone in it
                        if (v1type1 && degree[v1] == 1){
                            if(firstedge[v1]->end == v){
                                if (firstedge[v1]->invers->next->end != w){
                                    return FALSE;
                                }
                            } else if (firstedge[v1]->end == w){
                                if (firstedge[v1]->invers->next->end != v){
                                    return FALSE;
                                }
                            } else {
                                return FALSE;
                            }
                            break;
                        } else {
                            return FALSE;
                        }
                    }
                    break;
                case 3: //v1 must be in one of the faces alone
                    if (degree[v1] == 1){
                        if(firstedge[v1]->end == v){
                            if (firstedge[v1]->invers->next->end != w){
                                return FALSE;
                            }
                        } else if (firstedge[v1]->end == w){
                            if (firstedge[v1]->invers->next->end != v){
                                return FALSE;
                            }
                        } else {
                            return FALSE;
                        }
                        break;
                    } else {
                        return FALSE;
                    }
                    break;
                default:
                    return FALSE; // more than 3 edges between v and w: no choice of vi is possible
            }
        }
    }
    return TRUE;
}

/**************************************************************************/
static int
mark_subgraph_2() {
    int marked = FALSE;
    EDGE* edge;
    EDGE* edge2;
    RESETMARKS;
    RESETMARKS_V;
    for (int i = 0; i < nv; ++i) {region[i] = -1;}

    //mark all parallel edges
    for (int v = 0; v < nv; ++v) {
        if (class[v] == 0) {
            edge = firstedge[v];
            //For each edge we find the first edge in the same cyclic next-orbit that has the same end. If it is not edge, then we have found a 2-cycle
            do {
                edge2 = edge->next; //edge->next is not a candidate, then there is no edge between edge and edge2
                while (edge->end != edge2->end) {
                    edge2 = edge2->next;
                }
                if (edge2 != edge) {
                    //parallel edges found
                    edge->index = 1; //indicates that we have not yet handled the face containing this edge
                    edge->invers->index = 1;
                    MARK(edge);
                    MARK(edge->invers);
                    MARK(edge2);
                    MARK(edge2->invers);
                    MARK_V(edge->start);
                    MARK_V(edge->end);
                    region[edge->start] = region[edge->end] = 0;

                    marked = TRUE;
                }
                edge = edge->next;
            } while (edge != firstedge[v]);
        }
    }
    return marked;
}


/**************************************************************************/
static void
assign_c2_regions() {
    //assumes that all parallel edges have been marked
    //returns the number of regions that need to have a vi. region[v] is the number of the region that contains
    //if no region has been assigned yet, it is -1
    //if a region does not need a vi, it is 0
    EDGE* edge;
    EDGE* edge2;

    if (!mark_subgraph_2()){
        nregions = 0;
        return;
    }

    int nextregion = 1;
    int succes;
    for (int v = 0; v < nv; ++v) {
        edge = firstedge[v];
        //For each edge fill its corresponding face if necesary
        do {
            if(ISMARKED(edge) && edge->index){
                if (marked_facesize(edge) == 2){
                    succes = assign_region(edge, nextregion);
                    if (succes == nextregion){
                        nextregion++;
                    } else if (succes > 0){ //this region should be made 0
                        for (int w = 0; w < nv; ++w) {
                            if (region[w] == succes || region[w] == nextregion){
                                region[w] = 0;
                            } else if (region[w] > succes){
                                region[w] = region[w] - 1;
                            }
                        }
                        nextregion--;
                    }
                } else {
                    succes = assign_region(edge, 0);
                    if (succes != 0){//this region should be made 0
                        nextregion--;
                        for (int w = 0; w < nv; ++w) {
                            if (region[w] == succes){
                                region[w] = 0;
                            } else if (region[w] > succes){
                                region[w] = region[w] - 1;
                            }
                        }
                    }
                }
                unindex_markedface(edge);
            }
            edge = edge->next;
        } while (edge != firstedge[v]);
    }
    nregions = nextregion - 1;
}

/**************************************************************************/
static int is_c3(uint16_t v0, uint16_t v1, uint16_t v2, int v1type1){
    int x, y;
    int intersection[nv];
    int intersectioncount;

    EDGE *edge, *edgestart;

    if (nv == 3){//P3 is a special case when it comes to cycles, so it is treated separately
        return TRUE;
    }

    for (int i = nregions; i > 0; --i) {
        //every region must be represented
        if (!(v1type1 && i == region[v1]) && !(region[v0] == i || region[v2] == i)) {
            return FALSE;
        }
    }

    for (int v = 0; v < nv; ++v)
        if (class[v] == 0)
            for (int w = v; w < nv; ++w)
                if (class[w] == 0) {
                    //v and w loop over all pairs of vertices of class 0. v == w is allowed
                    intersectioncount = 0;
                    for (x = 0; x < nv; ++x) {
                        if (adjacency[v][x] && adjacency[w][x]) {
                            intersection[intersectioncount++] = x;
                        }
                    }
                    if (intersectioncount > 0) {
                        //the intersection is not empty, there may be non-empty 4-cycles or 2-cycles

                        //we loop over all pairs of vertices in the intersection of the neighbours of v and w, x == y is allowed
                        for (int intersect1 = 0; intersect1 < intersectioncount; ++intersect1) {
                            x = intersection[intersect1];
                            for (int intersect2 = intersect1; intersect2 < intersectioncount; ++intersect2) {
                                y = intersection[intersect2];

                                //the subgraph induced by v, w, x and y is now marked
                                RESETMARKS;
                                RESETMARKS_V;
                                MARK_V(v);
                                MARK_V(w);
                                MARK_V(x);
                                MARK_V(y);
                                mark_edges_to_marked_vertices(x);
                                mark_edges_to_marked_vertices(y);

                                if (v == w && x == y && v1type1) { //then v1 has degree 1 and it can be ignored. It is as if its 2-cycle is empty
                                    MARK_V(v1);
                                }

                                //find 2- and 4-cycles
                                edgestart = firstedge[v];
                                edgestart = next_marked_edge(edgestart);
                                edge = edgestart;
                                do {
                                    if (next_marked_edge(next_marked_edge(edge->invers)->invers) == edge) {
                                        //a face of size 2
                                        if (v == w && x == y) { //only check inside once, when we are looking at the edges between v==w and x==y
                                            if (!(vi_inside_or_empty(edge, v0, nv, v2))) {
                                                return FALSE;
                                            }
                                        }
                                    } else {
                                        //a face of size 4
                                        if (!(vi_inside_or_empty(edge, v0, v1, v2))) {
                                            return FALSE;
                                        }
                                    }
                                    edge = next_marked_edge(edge);
                                } while (edge != edgestart);
                            }
                        }
                    }
    }
    return TRUE;
}

/**************************************************************************/
static int
is_c2(uint16_t v0, uint16_t v1, uint16_t v2) {
    if (nregions > 3){
        return FALSE;
    }

    for (int i = nregions; i > 0; --i) {
        //every region must be represented
        if (!(region[v0] == i || region[v1] == i || region[v2] == i)) {
            return FALSE;
        }
    }
    return TRUE;
}

/**************************************************************************/
static int
connectivity(uint16_t v0, uint16_t v1, uint16_t v2, int v1type1) {
    if (!is_c2(v0,v1,v2)){
        return 1;
    } else if ((min_connectivity == 2 && !sswitch) || !maybe_c3(v0,v1,v2, v1type1) || !is_c3(v0,v1,v2, v1type1)){
        return 2;
    }
    return 3;
}


/**************************************************************************/
static void
write_lopsp(uint16_t v0, uint16_t v1, uint16_t v2, uint16_t v0type, uint16_t v1type1, int mirror){
    /* writes the lopsp-code of the current lopsp-operation. Assumes edgecode has been computed */

    if (uswitch){
        return;
    }
    uint16_t v0_print = v0 + 1;
    uint16_t v1_print = v1 + 1;
    uint16_t v2_print = v2 + 1; //in lopsp-format vertices count from 1, not 0

    fwrite(&nv, sizeof(uint16_t), 1, outfile);
    fwrite(&v0_print, sizeof(uint16_t), 1, outfile);
    fwrite(&v1_print, sizeof(uint16_t), 1, outfile);
    fwrite(&v2_print, sizeof(uint16_t), 1, outfile);
    fwrite(&v0type, sizeof(uint16_t), 1, outfile);
    fwrite(&v1type1, sizeof(uint16_t), 1, outfile);

    fwrite(mirror ? edgecode_mirror : edgecode,sizeof(uint16_t),edgecode_length,outfile);
}


/**************************************************************************/
static void
compute_code(){
    /* computes the edgecode part of the lopsp-encoding for lopsp-operations associated with the current quadrangulation.*/
    register EDGE *edge;
    int i,j;
    uint16_t* code = edgecode; //points to the current place in edgecode
    uint16_t edgenumber;
    RESETMARKS

    edgecode_length = nv + ne;//the number of zeroes plus the number of directed edges
    edgenumber = 1;

    for (i = 0; i < nv; i++){
        for (j = 0, edge = firstedge[i]; j < degree[i]; j++, edge = edge->next){
            if (ISMARKED(edge)) { *code = edge->index;}
            else { //a marked edge already got its number
                MARK(edge->invers);
                edge->invers->index = edgenumber;
                *code = edgenumber;
                edgenumber++;
            }
            code++;
        }
        *code = 0;  code++;
    }
}

/**************************************************************************/
static void
compute_code_mirror(){
    /* computes the edgecode part of the lopsp-encoding for the mirror image of the lopsp-operations associated
     * with the current quadrangulation.*/
    register EDGE *edge;
    int i,j;
    uint16_t* code = edgecode_mirror; //points to the current place in edgecode_mirror
    uint16_t edgenumber;
    RESETMARKS

    edgenumber = 1;

    for (i = 0; i < nv; i++){
        for (j = 0, edge = firstedge[i]; j < degree[i]; j++, edge = edge->prev){
            if (ISMARKED(edge)) { *code = edge->index;}
            else { //a marked edge already got its number
                MARK(edge->invers);
                edge->invers->index = edgenumber;
                *code = edgenumber;
                edgenumber++;
            }
            code++;
        }
        *code = 0;  code++;
    }
}

/**************************************************************************/
static int
are_canonical_vi(int number_of_automorphisms, int v0, int v1, int v2, const int* certain_vi, int certain_vicount, int v1type1){
    /* Returns 0 if this choice of vi is not valid. A valid choice is canonical and has the vi in certain_vi
     * Returns 1 if switching v0 and v2 gives an isomorphic lopsp-operation, 2 if it does not. */
    //the canonical form has the lexicographically highest numbers

    if (degree[v1] == 1 && !v1type1 && min_connectivity == 3){
        return FALSE;
    }

    //check if all certain_vi are chosen
    for (int i = 0; i < certain_vicount; ++i) {
        if (certain_vi[i] != v0 && certain_vi[i] != v1 && certain_vi[i] != v2){
            return FALSE;
        }
    }

    //check if canonical
    if (number_of_automorphisms == 1){
        return 2;
    }

    int other_v0, other_v1, other_v2;
    int to_return = 2;

    for (int automorphism = 0; automorphism < number_of_automorphisms; ++automorphism) {
        other_v0 = numbering[automorphism][vindex[v0]]->start;
        other_v1 = numbering[automorphism][vindex[v1]]->start;
        other_v2 = numbering[automorphism][vindex[v2]]->start;
        if (other_v0 < other_v2) {
            if (other_v0 > v0) {
                return FALSE;
            }
            if (other_v0 == v0) {
                if (other_v1 > v1) {
                    return FALSE;
                }
                if (other_v1 == v1) {
                    if (other_v2 > v2) {
                        return FALSE;
                    }
                }
            }
        } else if (other_v0 > other_v2) {
            if (other_v2 > v0) {
                return FALSE;
            }
            if (other_v2 == v0) {
                if (other_v1 > v1) {
                    return FALSE;
                }
                if (other_v1 == v1) {
                    if (other_v0 > v2) {
                        return FALSE;
                    }
                }
            }
        }
        if(other_v0 == v2 && other_v1 == v1 && other_v2 == v0){//there is an automorphism that switches v0 and v2
            to_return = 1;
        }
    }

    return to_return; //no bigger found
}
/**************************************************************************/
static int
mirror_different(int v0, int v1, int v2, int nbop, int nbtot){
    //Check for all orientation-reversing automorphisms if they map this lopsp-operation to its mirror image
    //if so, then the mirror image must not be written and the lopsp-operation is an lsp-operation
    if (nbop == nbtot || nbop == 0){
        return TRUE;
    }

    for (int automorphism = nbop; automorphism < nbtot; ++automorphism) {
        if (numbering[automorphism][vindex[v0]]->start == v0
         && numbering[automorphism][vindex[v1]]->start == v1
         && numbering[automorphism][vindex[v2]]->start == v2) {
            //there is an equivalent mirror copy
            return FALSE;
        }
    }

    return TRUE;
}

/**************************************************************************/
static int
count_lopsps_v1type02(int nbop, int nbtot, const int* possible_vi, const int* certain_vi, int certain_vicount){
    int v0,v1,v2; //current choices.
    int connectivity_quad, lsp, switchv02;
    int inf_factor = (2 * nv) - 4;
    int count = 0;

    if (lswitch && (nbop == nbtot || nbop == 0)){ //then no lsp-operations
        return 0;
    }
    if (nbtot == 1 && (min_connectivity == 1 || (min_connectivity == 2 && nregions == 0)) && uswitch && !sswitch){
        nout_op[inf_factor][1] += 4 * nv * (nv - 1) * (nv - 2);
        nout[inf_factor][1] += 2 * nv * (nv - 1) * (nv - 2);
        return oswitch ? 4 * nv * (nv - 1) * (nv - 2) : 2 * nv * (nv - 1) * (nv - 2);
    }

    for (v0 = 0; v0 < nv; v0++)
        if (possible_vi[v0])
            for (v1 = 0; v1 < nv; v1++)
                if (v1 != v0 && possible_vi[v1])
                    for (v2 = v0 + 1; v2 < nv; v2++)
                        if (v2 != v1 && v2 != v0 && possible_vi[v2]) //here we loop over all possible choices of v0,v1,v2 with v0 < v2
                            if ((switchv02 = are_canonical_vi(nbtot, v0, v1,
                                                              v2, certain_vi, certain_vicount, FALSE))) { //determine if canonical. Every lopsp has exactly one canonical choice
                                if (sswitch || min_connectivity > 1){
                                    connectivity_quad = connectivity(v0, v1, v2, FALSE);
                                } else {
                                    connectivity_quad = 1;
                                }
                                if (connectivity_quad >= min_connectivity) {
                                    lsp = !mirror_different(v0, v1, v2, nbop, nbtot);
                                    if (!lsp && !lswitch) {
                                        // not lsp, so not printed if lswitch
                                        //print mirror images
                                        if (oswitch) {
                                            write_lopsp(v0, v1, v2, 0, 0, 1);
                                            write_lopsp(v0, v1, v2, 2, 0, 1);
                                            if (switchv02 == 2){
                                                write_lopsp(v2, v1, v0, 0, 0, 1);
                                                write_lopsp(v2, v1, v0, 2, 0, 1);
                                            }
                                        }
                                        nout_op[inf_factor][connectivity_quad] += switchv02 * 2;
                                        count += oswitch ? switchv02 * 2 : 0;
                                    }
                                    if ((lswitch && lsp) || !lswitch) {
                                        //if lswitch, only print if lsp, else always print
                                        write_lopsp(v0, v1, v2, 0, 0, 0);
                                        write_lopsp(v0, v1, v2, 2, 0, 0);
                                        if (switchv02 == 2){
                                            write_lopsp(v2, v1, v0, 0, 0, 0);
                                            write_lopsp(v2, v1, v0, 2, 0, 0);
                                        }
                                        nout_op[inf_factor][connectivity_quad] += switchv02 * 2;
                                        nout[inf_factor][connectivity_quad] += switchv02 * 2;
                                        count += switchv02 * 2;
                                    }
                                }
                            }
    return count;
}

/**************************************************************************/
static int
count_lopsps_v1type1(int nbop, int nbtot, const int* possible_vi, const int* certain_vi, int certain_vicount){
    int v0,v1,v2; //current choices.
    int connectivity_quad, lsp, switchv02;
    int inf_factor = (2 * nv) - 5;
    int count = 0;

    if (lswitch && (nbop == nbtot || nbop == 0)){ //then no lsp-operations
        return 0;
    }
    if (nbtot == 1){
        if (min_connectivity == 1 && uswitch && !sswitch){
            int deg1count = 0;
            for (int i = 0; i < nv; ++i) {
                if (degree[i] == 1){
                    deg1count++;
                }
            }
            nout_op[inf_factor][1] += 2 * deg1count * (nv - 1) * (nv - 2);
            nout[inf_factor][1] += deg1count * (nv - 1) * (nv - 2);
            return oswitch ? 2 * deg1count * (nv - 1) * (nv - 2) : deg1count * (nv - 1) * (nv - 2);
        }
    }

    for (v1 = 0; v1 < nv; v1++)
        if (degree[v1] == 1 && possible_vi[v1]) {
            for (v0 = 0; v0 < nv; v0++)
                if (v1 != v0 && possible_vi[v0])
                    for (v2 = v0 + 1; v2 < nv; v2++)
                        if (v2 != v1 && v2 != v0 &&
                            possible_vi[v2]) //here we loop over all possible choices of v0,v1,v2
                            if ((switchv02 = are_canonical_vi(nbtot, v0, v1, v2, certain_vi,
                                                              certain_vicount, TRUE))) { //determine if canonical
                                if (sswitch || min_connectivity > 1) {
                                    connectivity_quad = connectivity(v0, v1, v2, TRUE);
                                } else {
                                    connectivity_quad = 1;
                                }
                                if (connectivity_quad >= min_connectivity) {
                                    lsp = !mirror_different(v0, v1, v2, nbop, nbtot);
                                    if (!lsp && !lswitch) {
                                        // not lsp, so not printed if lswitch
                                        //print mirror imagesw
                                        if (oswitch) {
                                            write_lopsp(v0, v1, v2,
                                                        (class[firstedge[v1]->end] == class[v0]) ? 0 : 2, 1, 1);
                                            if (switchv02 == 2) {
                                                write_lopsp(v2, v1, v0,
                                                            (class[firstedge[v1]->end] == class[v2]) ? 0 : 2, 1, 1);
                                            }
                                        }
                                        nout_op[inf_factor][connectivity_quad] += switchv02;
                                        count += oswitch ? switchv02 : 0;
                                    }
                                    if ((lswitch && lsp) || !lswitch) {
                                        //if lswitch, only print if lsp, else always print
                                        write_lopsp(v0, v1, v2, (class[firstedge[v1]->end] == class[v0]) ? 0 : 2, 1,
                                                    0);
                                        if (switchv02 == 2) {
                                            write_lopsp(v2, v1, v0,
                                                        (class[firstedge[v1]->end] == class[v2]) ? 0 : 2, 1, 0);
                                        }
                                        nout_op[inf_factor][connectivity_quad] += switchv02;
                                        nout[inf_factor][connectivity_quad] += switchv02;
                                        count += switchv02;
                                    }
                                }
                            }
        }
    return count;
}

/**************************************************************************/
static int
c2_candidate(){
    int deg1count = 0;
    for (int v = 0; v < nv; ++v) {
        if (degree[v] == 1) {
            deg1count++;
            if (deg1count > 3) {
                return FALSE;
            }
        }
        for (int w = v + 1; w < nv; ++w) {
            if (adjacency[v][w] > 3){
                return FALSE;
            }
        }
    }
    return TRUE;
}

/**************************************************************************/

static void
got_quadrangulation(int nbtot, int nbop)

/* This is called when a quadrangulation of maxnv vertices is found.  The main
   purpose is to write the graph and to collect some stats. */
{
    int possible_vi[nv];
    int certain_vi[3];
    int deg1count = 0;
    int deg2count = 0;
    int lopsp_count = 0;

    total_quads_generated++;

    //possible_vi[v] will be TRUE if vertex v may be one of the vi
    for (int v = 0; v < nv; ++v) {
        possible_vi[v] = TRUE;
    }

    if (min_connectivity >= 2) {
        if (min_connectivity == 2) {
            for (int v = 0; v < nv; ++v) {
                if (degree[v] == 1) {
                    deg1count++;
                    if (deg1count > 3) {
                        return;
                    }
                    certain_vi[deg1count - 1] = v;
                }
            }
        } else if (min_connectivity == 3 && nv > 4) { //for nv = 4 not every vertex of degree 2 must be a vi
            for (int v = 0; v < nv; ++v) {
                if (degree[v] < 3) {
                    if (degree[v] == 1) {
                        deg1count++;
                    } else if (degree[v] == 2) {
                        deg2count++;
                    }
                    if (deg1count + deg2count > 3) {
                        return;
                    }
                    certain_vi[deg1count + deg2count - 1] = v;
                }
            }
        }
        for (int v = 0; v < nv; ++v) {
            for (int w = v + 1; w < nv; ++w) {
                switch (adjacency[v][w]) {
                    case 0:
                    case 1:
                    case 2: //not both v and w may be vi, but that is difficult to describe
                        break;
                    case 3:
                        possible_vi[v] = possible_vi[w] = FALSE;
                        break;
                    default:
                        return; // more than 3 edges between v and w: no choice of vi is possible
                }
            }
        }
    }

    compute_code(); //computes the edgecode. Will be written for every lopsp-operation associated with the quadrangulation
    if (oswitch){
        compute_code_mirror();
    }

    if (nbtot != 1){
        determine_edge_numbers();
    }
    if (sswitch || min_connectivity > 1) {
        if (sswitch && min_connectivity != 3){
            if (c2_candidate()){
                assign_c2_regions();
            } else {
                nregions = 4;//too high, will not be considered for c2
            }
        } else {
            assign_c2_regions();
            if (nregions > 3){
                return;
            }
        }
    }

    if (!aswitch){
        if (factor % 2 == 0){
            lopsp_count += count_lopsps_v1type02(nbop, nbtot, possible_vi, certain_vi, deg1count + deg2count);
        } else {
            lopsp_count += count_lopsps_v1type1(nbop, nbtot, possible_vi, certain_vi, deg1count + deg2count);
        }
    } else {
        lopsp_count += count_lopsps_v1type1(nbop, nbtot, possible_vi, certain_vi, deg1count + deg2count);
        if ((maxnv != nv || factor % 2 == 0) && !(min_connectivity == 3 && nregions > 2)){
            lopsp_count += count_lopsps_v1type02(nbop, nbtot, possible_vi, certain_vi, deg1count + deg2count);
        }
    }

    if (lopsp_count > 0){
        totalquads++;
    }
}


/****************************************************************************/
static void
scanmultiquadrangulations(int nbtot, int nbop)

/* The main node of the recursion for general quadrangulations.
   As this procedure is entered, nv,ne,degree etc are set for some graph,
   and nbtot/nbop are the values returned by canon() for that graph.
*/

{
    EDGE *firstedge_save[MAXN];
    EDGE *extP0[MAXE],*extP1[MAXE],*extD1[MAXE];
    EDGE *good_or[MAXE],*good_mir[MAXE];
    int ngood_or,ngood_mir,ngood_ref,ngood_mir_ref;
    int nextP0,nextP1,nextD1,i;
    int xnbtot,xnbop;
    EDGE *rededge;

    if (nv == maxnv)
    {
        got_quadrangulation(nbtot, nbop);  /* Third arg is connectivity but
                                       we don't know it. */
        return;
    } else if (aswitch){
        got_quadrangulation(nbtot, nbop);
    }

    if (nv == splitlevel)
    {
        if (splitcount-- != 0) return;
        splitcount = mod - 1;

        for (i = 0; i < nv; ++i) firstedge_save[i] = firstedge[i];
    }

    //finds all possible extensions
    find_extensions_multiquad(nbtot,nbop,extD1,&nextD1,extP0,&nextP0,extP1,&nextP1);

    for (i = 0; i < nextD1; ++i)//try all D1 extensions
    {
        rededge = extend_quadr_D1(extD1[i]);
        {
            quadr_D1_legal_all(rededge, good_or, &ngood_or, &ngood_ref, good_mir,
                               &ngood_mir, &ngood_mir_ref);
            if (ngood_ref + ngood_mir_ref > 0) {
                if (nv == maxnv && !needgroup && ngood_or == ngood_ref
                    && ngood_mir == ngood_mir_ref)
                    got_quadrangulation(1, 1);
                else if (ngood_or + ngood_mir == 1)
                    scanmultiquadrangulations(1, 1);
                else if (canon_edge_oriented(good_or, ngood_or, ngood_ref,
                                             good_mir, ngood_mir, ngood_mir_ref,
                                             degree, numbering, &xnbtot, &xnbop)) {
                    scanmultiquadrangulations(xnbtot, xnbop);
                }
            }
        }
        reduce_quadr_D1(rededge);
    }

    for (i = 0; i < nextP0; ++i)//then try all P0 extensions
    {
        rededge = extend_quadr_P0(extP0[i]);
        {
            if (!has_quadr_D1()) {
                multiquadr_P0_legal_all(rededge, good_or, &ngood_or, &ngood_ref, good_mir,
                                        &ngood_mir, &ngood_mir_ref);
                if (ngood_ref + ngood_mir_ref > 0) {
                    if (nv == maxnv && !needgroup && ngood_or == ngood_ref
                        && ngood_mir == ngood_mir_ref)
                        got_quadrangulation(1, 1);
                    else if (ngood_or + ngood_mir == 1)
                        scanmultiquadrangulations(1, 1);
                    else if (canon_edge_oriented(good_or, ngood_or, ngood_ref,
                                                 good_mir, ngood_mir, ngood_mir_ref,
                                                 degree, numbering, &xnbtot, &xnbop)) {
                        scanmultiquadrangulations(xnbtot, xnbop);
                    }
                }
            }
        }
        reduce_quadr_P0(rededge);
    }

    for (i = 0; i < nextP1; ++i)//finally try all P1 extensions
    {
        rededge = extend_quadr_P1(extP1[i]);
        {
            if (!has_quadr_D1() && !has_quadr_P0_multi())
            {
                multiquadr_P1_legal_all(rededge,good_or,&ngood_or,&ngood_ref,
                                        good_mir,&ngood_mir,&ngood_mir_ref);
                if (ngood_ref+ngood_mir_ref > 0
                    && canon_edge_oriented(good_or,ngood_or,ngood_ref,
                                           good_mir,ngood_mir,ngood_mir_ref,
                                           degree,numbering,&xnbtot,&xnbop)) {
                    scanmultiquadrangulations(xnbtot, xnbop);
                }
            }
        }
        reduce_quadr_P1(rededge);
    }

    if (nv == splitlevel)
        for (i = 0; i < nv; ++i) firstedge[i] = firstedge_save[i];
}



/****************************************************************************/

static void
open_output_file(void)

/* Open the output file, and write a header if one is called for.
   All the needed information is in global vars.
*/
{
    msgfile = stderr;
    outfilename = "stdout";
    outfile = stdout;

    if(!uswitch) { //if output is required
        fputs(">>lopsp<<", outfile);
    }
}

/****************************************************************************/

static void
initialize_splitting(int minlevel, int hint, int maxlevel)

/* Set splitlevel and splitcount.  minlevel and maxlevel are bounds
   on its value.  It must be that both minlevel and maxlevel are at
   least equal to the largest starting order (nv for external calls
   to scansimple() or similar routines), and at most equal to the
   smallest parent of an internal-output graph (call from scansimple()
   or similar to got_one() or similar).

   hint is a desirable value, which can be anything as the actual
   value used is forced between minlevel and maxlevel.  For plugins,
   the value of splithint is used instead if it is >= 0.

   In case there is no way to use splitting within those limits,
   it is turned off by setting splitlevel=0.  In that case only
   subcase 0 should produce output.

   Splitting occurs at the first level where nv >= splitlevel.
   If an operation can add k vertices, it must be that
   splitlevel <= maxnv - k.
*/
{
    splitlevel = hint ;

    if (splitlevel > maxlevel) splitlevel = maxlevel;

    if (splitlevel < minlevel && splitlevel > 0)
    {
        if (minlevel <= maxlevel) splitlevel = minlevel;
        else                      splitlevel = 0;
    }
    if (mod == 1) splitlevel = 0;

    splitcount = res;
}
/****************************************************************************/

static void
multiquadrangulation_dispatch(void)

/* General quadrangulations. */
{
    int startingsize,nbtot,nbop,hint;

    startingsize = 3;

    open_output_file();

    needgroup = 1;

    hint = (maxnv < 16 ? 12 : 13);
    initialize_splitting(startingsize,hint,maxnv-4);
    if (splitlevel == 0 && res > 0) return;
    initialize_multiquadrangulations();
    canon(degree,numbering,&nbtot,&nbop);        // makes initial group
    scanmultiquadrangulations(nbtot,nbop);
}

/****************************************************************************/

static int
getswitchvalue(const char *arg, int pj)

/* Find integer value for switch.
   arg is a pointer to a command-line argument.
   pj is an index into arg, which is updated.
   The value of the switch is the function return value.
   For example, if arg="-xyz1432q" and *pj=3 (pointing at "z"),
       the value 1432 is returned and *pj=7 (pointing at "2").
   An absent value is equivalent to 0.
*/

{
    int j,ans;

    ans = 0;
    for (j = pj; arg[j+1] >= '0' && arg[j+1] <= '9'; ++j)
        ans = ans * 10 + (arg[j+1] - '0');

    return ans;
}

/****************************************************************************/

static int
getswitchvalue_and_move(const char *arg, int* pj)

/* Find integer value for switch.
   arg is a pointer to a command-line argument.
   pj is an index into arg, which is updated.
   The value of the switch is the function return value.
   For example, if arg="-xyz1432q" and *pj=3 (pointing at "z"),
       the value 1432 is returned and *pj=7 (pointing at "2").
   An absent value is equivalent to 0.
*/

{
    int ans = 0;

    while (arg[(*pj)+1] >= '0' && arg[(*pj)+1] <= '9'){
        ans = ans * 10 + (arg[(*pj)+1] - '0');
        (*pj)++;
    }
    return ans;
}

/****************************************************************************/

static void
decode_command_line(int argc, char *argv[])

/* Decode the command line, setting the global variables which
   give the switch values.  Some basic checking is done too, but
   the most detailed checking is done later.  If an error is
   found, this procedure never returns.
*/

{
    int index;
    lswitch = uswitch = oswitch = sswitch = FALSE;

    res = 0; mod = 1;
    min_connectivity = 1;

    int option;
    while ((option = getopt(argc, argv, SWITCHES)) != -1) {
        switch (option) {
            case 'l': //only lsp-operations
                lswitch = TRUE;
                break;
            case 'a': //all operations with inflation factor at most nv
                aswitch = TRUE;
                break;
            case 'h': //help
                print_help();
                exit(0);
                break;
            case 'u': //do not write the operations
                uswitch = TRUE;
                break;
            case 'o': //only consider orientation preserving automorphisms
                oswitch = TRUE;
                break;
            case 'm': //res/mod pair
                index = -1;
                res = getswitchvalue_and_move(optarg, &index);
                if (optarg[index + 1] != '/'){
                    fprintf(stderr, "The option -m needs a value of the form res/mod, where res and mod are non-negative integers such that res < mod.\n");
                    exit(1);
                }
                index++;
                mod = getswitchvalue(optarg, index);

                if (res >= mod){
                    fprintf(stderr, "The option -m needs a value of the form res/mod, where res and mod are non-negative integers such that res < mod.\n");
                    exit(1);
                }
                break;
            case 'c': // set minimum connectivity to 2 or 3
                if (optarg[0] == '1' || optarg[0] == '2' || optarg[0] == '3'){
                    min_connectivity = optarg[0] - '0';
                } else {
                    fprintf(stderr, "The option -c must be followed by 1 (default), 2 or 3. Other values are not supported.\n");
                    exit(1);
                }
                if (optarg[1] != '\0'){
                    fprintf(stderr, "The option -c must be followed by 1 (default), 2 or 3. Other values are not supported.\n");
                    exit(1);
                }
                break;
            case 's': //write various stats
                sswitch = TRUE;
                break;
            case ':':
                fprintf(stderr, "option needs a value.\n");
                exit(1);
                break;
            default:
                fprintf(stderr, "unrecognized option %c\n", option);
                exit(1);
        }
    }

    if (optind >= argc){
        fprintf(stderr, "This program always takes a number as input. If you are using switches, check that they have a value if they require one.\n");
        exit(1);
    }

    if (aswitch && mod > 1){
        fprintf(stderr, "The options -a and -m cannot be used together.\n");
        exit(1);
    }

    factor = getswitchvalue(argv[optind], -1);
    /* If v1 is not of type 1, The inflation factor is the number of edges of the quadrangulation: 2*maxnv - 4
    If v1 is of type 1. The inflation factor is the number of edges of the quadrangulation minus one: 2*maxnv - 5
    */
    maxnv = (factor % 2 == 0) ? (factor + 4)/ 2 : (factor + 5)/ 2;

    outfilename = "stdout";
}

/****************************************************************************/
int ndigits(bigint number){
    bigint bignumber = number;
    int count = 1;
    while (bignumber >= 10){
        count++;
        bignumber = bignumber / 10;
    }
    return count;
}

/****************************************************************************/
int
main(int argc, char *argv[])
{
    int i;
    decode_command_line(argc, argv);

    struct tms timestruct0,timestruct1;

    times(&timestruct0);

    if (aswitch){
        for(int f = 0; f < factor; f++){
            for (i = 0; i < 6; ++i){
                nout[f][i] = 0;
                nout_op[f][i] = 0;
            }
        }
    }
    for (i = 0; i < 6; ++i) nout[factor][i] = nout_op[factor][i] = 0;
    totalquads = 0;

    multiquadrangulation_dispatch();

    times(&timestruct1);

    totalout = totalout_op = 0;
    for (i = 0; i < 6; ++i)
    {
        totalout += nout[factor][i];
        totalout_op += nout_op[factor][i];
    }
    if (aswitch){
        for(int f = 0; f < factor; f++){
            for (i = 0; i < 6; ++i){
                totalout += nout[f][i];
                totalout_op += nout_op[f][i];
            }
        }
    }
    if(mod > 1){
        fprintf(msgfile,"Res/mod splitting: part %d modulo %d.\n", res, mod);
    }

    PRINTBIG(msgfile,(oswitch ? totalout_op : totalout));
    if (min_connectivity > 1){
        fprintf(msgfile," %d-connected", min_connectivity);
    }
    fprintf(msgfile," Lopsp-operations of inflation factor %s%d", aswitch? "at most ": " ",factor);
    if (uswitch) fprintf(msgfile," counted");
    else         fprintf(msgfile," written to %s",outfilename);
    fprintf(msgfile,", originating from %llu quadrangulations. ", totalquads);
    fprintf(msgfile,"; cpu=%.2f sec\n",
            (double)(timestruct1.tms_utime+timestruct1.tms_stime
                     -timestruct0.tms_utime+timestruct0.tms_stime) / (double)CLK_TCK);

    if (oswitch){
        fprintf(msgfile,"They are non-isomorphic up to orientation-preserving automorphism\n");
        fprintf(msgfile,"%llu of them are non-isomorphic up to any automorphism\n", totalout);
    }
    if (lswitch){
        fprintf(msgfile,"All of them can be written as an lsp-operation\n");
    }


    if (sswitch){ //print stats
        int width = MAX(ndigits(totalout), 15);
        for (int f = aswitch ? 1 : factor; f <= factor; f++){
            fprintf(msgfile,"\n");

            fprintf(msgfile,"f = %*d", 2, f);

            fprintf(msgfile,"|");
            if (!lswitch) {
                for (int j = 0; j < width - 11 + 2; ++j) {
                    fprintf(msgfile, " ");
                }
                fprintf(msgfile, "total count  |");
                for (int j = 0; j < width - 8 + 2; ++j) {
                    fprintf(msgfile, " ");
                }
                fprintf(msgfile, "total OP  |");
            }
            for (int j = 0; j < width - 3 + 2; ++j) {
                fprintf(msgfile, " ");
            }

            fprintf(msgfile,"lsp\n");

            if (min_connectivity == 1) {
                fprintf(msgfile, "  c1  |");
                if (!lswitch) {
                    fprintf(msgfile, "  %*llu  |", width, nout[f][1] + nout[f][2] + nout[f][3]);
                    fprintf(msgfile, "  %*llu  |", width, nout_op[f][1] + nout_op[f][2] + nout_op[f][3]);
                }
                fprintf(msgfile, "  %*llu  \n", width, 2 * (nout[f][1] + nout[f][2] + nout[f][3]) - (nout_op[f][1] + nout_op[f][2] + nout_op[f][3]));
            }

            if (min_connectivity <= 2) {
                fprintf(msgfile, "  c2  |");
                if (!lswitch) {
                    fprintf(msgfile, "  %*llu  |", width, nout[f][2] + nout[f][3]);
                    fprintf(msgfile, "  %*llu  |", width, nout_op[f][2] + nout_op[f][3]);
                }
                fprintf(msgfile, "  %*llu  \n", width, 2 * (nout[f][2] + nout[f][3]) - nout_op[f][2] - nout_op[f][3]);
            }

            if (min_connectivity <= 3) {
                fprintf(msgfile, "  c3  |");
                if (!lswitch) {
                    fprintf(msgfile, "  %*llu  |", width, nout[f][3]);
                    fprintf(msgfile, "  %*llu  |", width, nout_op[f][3]);
                }
                fprintf(msgfile, "  %*llu  \n", width, 2 * nout[f][3] - nout_op[f][3]);
            }
        }

    }

    return 0;
}
