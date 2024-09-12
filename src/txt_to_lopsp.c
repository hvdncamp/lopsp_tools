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

/* Converts a readable .txt-representation of a lopsp-operation to the .lopsp format.

Input is taken from stdin and output is given through stdout.

**THE TXT FORMAT**
The expected .txt-format is as follows. First the number n of vertices of the lopsp-operation is given. This is
followed by one line for every vertex (numbered 1 through n). Such a line consists of the number of the vertex,
followed by brackets ().These brackets contain the type (0, 1 or 2) of the vertex and if the vertex is v0, v1
or v2 this is also indicated. Then follows a colon ':' and the neighbours of the vertex separated by spaces. The
neighbours are given in the order of the embedding. An example for the operation gyro:

7
1(0,v2):2 6
2(1):1 6 3 6
3(0):2 6 4 7 6
4(1):3 6 5 6
5(0,v0):4 6
6(2):5 4 3 2 1 2 3 7 4
7(1,v1):3 6

Note that lopsp-operations may have double edges, so the embedding is not uniquely defined for we list adjacent
vertices instead of incident edges of a vertex. However, if an input encodes a plane triangulation then this is
the only plane embedding of the graph encoded by that input. This implies that there is only one way a given input
can be interpreted as a lopsp-operation.

**THE LOPSP FORMAT**
The .lopsp-format is a memory-efficient format to encode lopsp-operations. As every type 1 vertex of a
lopsp-operation O that is different from v1 has degree 4, we can remove it from the operation without
losing information. If v1 has type 1 we also remove one of its two edges. The resulting labeled embedded
graph O' (which is a plane quadrangulation) is what we will encode. The format starts with a header >>lopsp<<.
Then follows the number of vertices in O' and three unsigned shorts containing the numbers of the vertices
v0, v1 and v2 in that order. After that there is another unsigned short, which is 0 if v0 is of type 0 and 2
if v0 is of type 2. Then another, which is 0 if v1 is not of type 1 and 1 if it is. The types of these vertices
determine the types of all other vertices in O'. Finally, for every vertex its incident edges (every edge has a
unique number) are given in the order of the embedding, followed by a 0. Every number here is an unsigned short
(2 bytes). An example of the .lopsp-format for gyro:

>>lopsp<<53510110234506062153040

*/
#include<stdint.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define INITIAL_NEIGHBOURS_SIZE 32
#define FREE_NEIGHBOURS(n) {free((n)->list);}
#define MAX__VERTICES 256 // The maximum number of vertices the lopsp-operation can have.
#define BUFFERLENGTH (MAX__VERTICES * 5)

typedef struct Neighbourstruct{
    int size;
    int count;
    uint16_t* list;
} Neighbours;

void init_neighbours(int initialsize, Neighbours* ns){ 
    ns->size = initialsize;
    ns->count = 0;
    ns->list = malloc(initialsize* sizeof(uint16_t));
}

void duplicate_empty_neighbours(Neighbours* original, Neighbours* new){
    new->size = original->size;
    new->count = original->count;
    new->list = (uint16_t*)calloc(original->size, sizeof(uint16_t));
}

void add_neighbour(Neighbours* neighbours, uint16_t v){
    if(neighbours->count == neighbours->size){
        neighbours->list = (uint16_t*) realloc(neighbours->list, neighbours->size * 2 * sizeof(uint16_t));
    } 
    neighbours->list[neighbours->count] = v;
    neighbours->count++; 
}

//Checks if no angles appear more than once, otherwise prints an error message.
void unique_angles(Neighbours* neighbours, int count){
    for(int k=1; k<=count; k++){
        for (int i = 0; i < neighbours[k].count -1; i++){
            for(int j=i+2; j<neighbours[k].count -1; j++){
                if(neighbours[k].list[i]==neighbours[k].list[j] && neighbours[k].list[i+1] == neighbours[k].list[j+1]){
                    fprintf(stderr, "Angle %hd, %hd, %hd occurs more than once.\n", neighbours[k].list[i], k, neighbours[k].list[i+1]);
                    exit(1);
                }
            }   
            //the list is cyclic so the last and first vertex in the list also form a pair
            if(neighbours[k].list[i]==neighbours[k].list[neighbours[k].count -1] && neighbours[k].list[i+1] == neighbours[k].list[0]){
                fprintf(stderr, "Angle %hd, %hd, %hd occurs more than once.\n", neighbours[k].list[i], k, neighbours[k].list[i+1]);
                exit(1);
            }
        }
    }   
}

//Checks if first and second appear consecutively in the given list of neighbours, returns 0 if the angle was not found, else returns the index where the angle was found + 1
int has_angle(Neighbours* neighbours, uint16_t first, uint16_t second){ 
    for (int i = 0; i < neighbours->count -1; i++){
        if (neighbours->list[i]==first){
            if (neighbours->list[i+1] == second){
                return i+1 ;
            }
        }
    }
    
    if(neighbours->list[neighbours->count-1] == first){
        if(neighbours->list[0]==second){
            return neighbours->count;
        }
    } 
    return 0;
}

//prints the edgecode of the predecoration associated with the lopsp-operation.
void print_edgecode(Neighbours* neighbours, uint16_t n, const uint16_t* vi, const char* type){
    uint16_t index;
    uint16_t y;
    uint16_t zero = 0;

    Neighbours* copy = malloc((n+1)* sizeof(Neighbours));
    for (int i = 1; i <= n; i++)
    {
        duplicate_empty_neighbours(neighbours+i, copy +i);
    }  

    uint16_t count=0;

    int edge_decided = 0; //is 0 if it is decided which neighbour of v1 will be printed, 1 otherwise or if type of v1 is not 1

    for (uint16_t x = 1; x <= n; x++){
        if((type[x] != '1') || (vi[1] == x)){
            for(int i=0; i < neighbours[x].count; i++){
                y = neighbours[x].list[i];
                if(copy[x].list[i]!=0){ // we have seen this edge before
                    fwrite(&(copy[x].list[i]),sizeof(uint16_t), 1, stdout);
                } else if((type[x] != '1' && type[y] != '1') || ((vi[1] == y || x == vi[1]) && type[vi[1]] == '1' && edge_decided == 0)){ //we have not seen this edge before, and it is in the subgraph
                    index = has_angle(neighbours +y, neighbours[x].list[(i+1)%neighbours[x].count], x) % neighbours[y].count;
                    count++;
                    fwrite(&count,sizeof(uint16_t), 1, stdout);
                    copy[y].list[index] = count; //mark the edge with its new name
                    if (vi[1] == y || x == vi[1]){
                        edge_decided = 1;
                    }
                }
            }
            fwrite(&zero, sizeof(uint16_t), 1, stdout);
        }
    }
    for (int i = 1; i <= n; i++){
        FREE_NEIGHBOURS(copy +i)
    }
    free(copy);
}

int main(int argc,char *argv[]) {
    uint16_t typev0;
    uint16_t v1type1;
    uint16_t n;//Will contain the number of vertices in the lopsp-operation
    uint16_t curr; //The vertex number of the line that is currently being read
    uint16_t temp;
    uint16_t newcount=0;//Counts the number of vertices in the subgraph
    int edge_count;//Counts the number of directed edges
    uint16_t vi[3] ={0};//At index i is the vertex vi
    uint16_t y; //used for checking if the graph is a triangulation
    uint16_t z;
    uint16_t count;//a counter

    uint16_t* newnames = malloc((MAX__VERTICES + 1) * sizeof(uint16_t));//At index i is the new number of vertex i in the smaller graph

    Neighbours* neighbours = malloc((MAX__VERTICES + 1) * sizeof(Neighbours));//At index i is the list of neighbours of vertex i
    char* type = malloc((MAX__VERTICES + 1) * sizeof(char));//At index i is the type of vertex i
    char* neighbourbuffer = malloc(BUFFERLENGTH * sizeof(char));//to read the neighbours
    char* tok;

    fputs(">>lopsp<<", stdout);

    while (!feof(stdin)){ //Read one lopsp-operation in every iteration
        if (1 != fscanf(stdin, "%hd", &n)){
            fprintf(stderr, "The description of every lopsp-operation should start with a number.\n");
            exit(1);
        }
        //INPUT---------------------------------------------------------------------
        count = getc(stdin);
        if(count!='\n'){
            fprintf(stderr, "There should be a new line directly after the number of vertices of a lopsp-operation\n");
            exit(1);
        }

        curr=1; 
        newcount=0;
        edge_count=0;

        for (int i = 0; i <= n; i++){
            init_neighbours(INITIAL_NEIGHBOURS_SIZE, neighbours +i);
            newnames[i]=0;
        }

        while(curr<=n){//Reads one line in every iteration, curr is the vertex that is described in the current line
            if(scanf("%hd", &temp)!=1){
                fprintf(stderr, "The line that should describe vertex %hd does not start with a number\n", curr);
                exit(1);
            }

            if(temp!=curr){
                fprintf(stderr, "The line that should describe vertex %hd starts with number %hd\n", curr, temp);
                exit(1);
            }

            if(getc(stdin)!='('){
                fprintf(stderr, "no ( where it was expected at vertex %hd\n", curr);
                exit(1);
            }
            type[curr] = getc(stdin);
            if (type[curr] != '0' && type[curr] != '1' && type[curr] !='2'){
                fprintf(stderr, "%c is not a valid type for vertex %hd. It should be 0, 1 or 2\n", type[curr], curr);
                exit(1);
            }

            switch (getc(stdin))
            {
            case ',': //curr is one of the vi
                if(getc(stdin)!='v'){
                    fprintf(stderr, "There is a , that is not followed by a v at vertex %hd\n", curr);
                    exit(1);
                } 
                switch (getc(stdin))
                {
                case '0'://v0
                    typev0 = type[curr] - '0';
                    if (typev0 == 1){
                        fprintf(stderr, "v0 cannot be of type 1\n");
                        exit(1);
                    }
                    vi[0]=curr;                
                    break;
                case '1'://v1
                    v1type1 = type[curr] =='1' ? 1 : 0;
                    vi[1]=curr;
                    break;
                case '2'://v2
                    if (type[curr] =='1'){
                        fprintf(stderr, "v2 cannot be of type 1\n");
                        exit(1);
                    }
                    vi[2]=curr;
                    break;
                default:
                    fprintf(stderr, "There is a v without 0, 1 or 2 at vertex %hd\n", curr);
                    exit(1);
                }
                newcount++;
                newnames[curr]=newcount;
                if(getc(stdin)!=')'){
                    fprintf(stderr, "No ) where it was expected at vertex %hd\n", curr);
                    exit(1);
                }
                break;
            case ')':// curr is not one of the vi
                if(type[curr]=='0'||type[curr]=='2'){
                    newcount++;
                    newnames[curr]=newcount;
                } else if(type[curr]!='1'){
                    fprintf(stderr, "%c is not a valid type for vertex %hd. It should be 0, 1 or 2\n", type[curr], curr);
                    exit(1);
                }
                break;
            default:
                fprintf(stderr, "No ) where it was expected at vertex %hd. Perhaps there is a ' ' right before the ), or a ',' is missing.\n", curr);
                exit(1);
            }

            if(getc(stdin)!=':'){
                fprintf(stderr, "No : where it was expected at vertex %hd\n", curr);
                exit(1);
            }

            //Reading neighbours
            if(fgets(neighbourbuffer, BUFFERLENGTH, stdin) == NULL){//The rest of the line as a string
                fprintf(stderr, "The neighbours cannot be read. Are you sure your file is encoded correctly?\n");
                exit(1);
            }
            tok = strtok(neighbourbuffer, " \t\r");
            do{
                if(sscanf(tok, "%hd", &temp)==1){ //A number is read, everything else is ignored
                    if(temp > n){
                        fprintf(stderr, "There are only %hd vertices. Vertex %hd does not exist\n", n, temp);
                        exit(1);
                    } else if (temp <= curr && type[curr]==type[temp]){
                        fprintf(stderr, "There cannot be edges between vertices of the same type. Vertices %hd and %hd have the same type.\n", curr, temp);
                        exit(1);
                    } 
                    add_neighbour(neighbours + curr, temp);
                }       
            } while ((tok = strtok(NULL, " \t\r"))); // Continue until the end of the line

            //Check degree conditions for lopsp-operations
            if(neighbours[curr].count==0){
                fprintf(stderr, "Vertex %hd has no neighbours\n", curr);
                exit(1);
            }

            if(type[curr]=='1'){
                if(vi[1]==curr && neighbours[curr].count!=2){
                    fprintf(stderr, "Vertex %hd is v1 and is of type 1 so it should have degree 2\n", curr);
                    exit(1);
                } else if (vi[1]!=curr && neighbours[curr].count!=4){
                    fprintf(stderr, "Vertex %hd is not v1 and is of type 1 so it should have degree 4\n", curr);
                    exit(1);
                }
            }

            edge_count += neighbours[curr].count;
            curr++;
        }

        if(getc(stdin)!='\n' && !feof(stdin)){
            fprintf(stderr, "The line after vertex %hd should be empty\n", curr-1);
            exit(1);
        }

        //CHECKS--------------------------------------------------------------------------
        unique_angles(neighbours, n); //Check if all angles are unique

        for (int i = 0; i < 3; i++){
            if(vi[i]<=0){
                fprintf(stderr, "v%d was not marked in the .txt file. Read the description of the format for details on how to indicate which vertices are v0, v1 and v2\n", i);
                exit(1);
            }
            if(vi[i]>n){
                fprintf(stderr, "v%d was given a number higher than the total number of vertices in the graph. This should not be possible\n", i);
                exit(1);
            }
        }      

        //Check if all faces are triangles
        for (uint16_t x = 1; x <= n; x++){
            for(int i=0; i<neighbours[x].count; i++){
                y = neighbours[x].list[i];
                z = neighbours[x].list[(i+1)%neighbours[x].count];
                if (!has_angle(neighbours + y, z, x) || !has_angle(neighbours + z, x, y)){
                    fprintf(stderr, "The face containing the angle %hd, %hd, %hd is not a triangle\n", y, x, z);
                    exit(1);
                }
            }
        }

        //Check planarity. Using 2|E|=3|F| (which is true because all faces are triangles) we get that the graph is plane if and only if 6|V| = 12 + 2|E|. As every edge has two directed edges 2|E|=edge_count
        if (6*n != 12 + edge_count){
            fprintf(stderr, "The graph is not plane, |V|=%hd, |E|=%hd, |F|=%hd\n", n, edge_count/2, edge_count/3);
            exit(1);
        }

        //OUTPUT----------------------------------------------------------------------------------------------

        // Number of vertices
        fwrite(&newcount, sizeof(uint16_t), 1, stdout); 

        // Vertices v0, v1 and v2
        for (int i = 0; i < 3; i++){
            fwrite(newnames +(vi[i]), sizeof(uint16_t), 1, stdout);
        }

        // Type of v0
        fwrite(&typev0, sizeof(uint16_t), 1, stdout);
        fwrite(&v1type1, sizeof(uint16_t), 1, stdout);

        //edgecode
        print_edgecode(neighbours, n, vi, type);  

        for (int i = 0; i < n+1; i++){
            FREE_NEIGHBOURS(neighbours +i)
        }
    }

    free(newnames);
    free(neighbours);
    free(type); 
    free(neighbourbuffer);

    return(0);
}
