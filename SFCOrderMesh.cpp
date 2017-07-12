
//Brief Explaination
//######

//This is the file for space sapce filling curve implementation of the waveguide two 
//dimensional mesh. This will be documented well for the reading purpose. 

//Input : This file takes unit_circle.txt as a input. 
//Ouput : Mesh file with SFC ordering. 
//######

//# Header files 
#include <iostream>
#include <stdio.h>
#include <vector>
#include <stdint.h>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>
#include <iterator>
#include <chrono>


// Intel TBB library 
#include <tbb/parallel_sort.h>
#include "tbb/task_scheduler_init.h"

// OpenMP
#include <omp.h>

// Global Constanst 
#define zero 0 

// Data Strucutres 
struct vertices {
    
	uint64_t index; 
	uint64_t sfc_index; 
	double x; 
	double y; 
        bool flag ;
};

struct element {
    	
    uint64_t v1 ; 
    uint64_t v2 ; 
    uint64_t v3 ; 
    double centroidx ;
    double centroidy ;
    long int index;
    uint64_t sfc_index; 
    
};

// Data Structure of Bounding Box
struct BoundingBox {
    
    double centre_x  ;
    double centre_y ;
    double box_max_x ;
    double box_max_y ;
    double box_min_x ;
    double box_min_y ;
};

// Compare Functions Delcarations and Definitions
bool compareVerticesX(const vertices &v1,const vertices &v2) {
    
    if (v1.x < v2.x) { return true ;}
    else return false ;
}

bool compareVerticesY(const vertices&v1, const vertices&v2){
    
    if (v1.y < v2.y) {return true ;}
    else return false ;
    
}

// Function Declarations
uint64_t Morton_Curve(element E,BoundingBox box,uint64_t N_level) ;
bool compareElement(const element &e1, const element &e2) ;
bool compareVertices(const vertices &v1, const vertices &v2) ;

int main(int argc, char **argv){
    
    // File pointers
    std::string _filepointers;
    std::string _filename = argv[1];
    
    // Bounding Box Variables
    BoundingBox box ;
    BoundingBox tempBox ; 
    std::vector<double> maxCords(2,0.0) ;
    //double _maxBoxDim ;
    
    //Global Variable to keep hold Vertices Number 
    int64_t _globalVerticesNumber = 0 ; 
    
    // Number of vertices and Elements 
    uint64_t _numVertices; 
    uint64_t _numElements;
    
    // Vectors of Elements and Vertices 
    std::vector<element> mesh ; 
    std::vector<vertices> vertices_array ;

    //Opening the file 
    std::fstream gomfile; 
    gomfile.open(_filename, std::ios::out | std::ios::in); 

    gomfile>>_numVertices ; 
    
    std::string tempstring ; 
    std::getline(gomfile, tempstring) ; 
    std::getline(gomfile, tempstring) ; 
    
    vertices_array.resize(_numVertices) ; 
    
    std::cout << " The number of vertices are := "<<_numVertices<< std::endl ; 
    
    // Reading the vertices 
    double x,y ; 
    uint64_t id ; 
    for (size_t ind = 0 ; ind < _numVertices ; ++ind ){
        
        gomfile >> id ; 
        gomfile >> x ; 
        gomfile >> y ; 
        
        vertices_array[ind].index = id + 1 ; 
        vertices_array[ind].sfc_index = 0 ; 
        vertices_array[ind].x = x ; 
        vertices_array[ind].y = y ;  
        vertices_array[ind].flag = false ; 
        
    }
    
    std::getline(gomfile, tempstring) ; 
    
    gomfile >> _numElements ; 
    std::getline(gomfile, tempstring) ; 
    std::getline(gomfile, tempstring) ; 
    
    mesh.resize(_numElements) ; 
    
    // Reading the elements data
    for (size_t ind = 0 ; ind < _numElements ; ++ind){
            
        mesh[ind].index = ind + 1  ; 
        mesh[ind].sfc_index = 0 ; 
        
        gomfile >> mesh[ind].v1 ; 
        gomfile >> mesh[ind].v2 ; 
        gomfile >> mesh[ind].v3 ; 
        
        // Centroid 
        mesh[ind].centroidx = 0.0 ; 
        mesh[ind].centroidy = 0.0 ; 
        
    }
    
    std::cout<<"The File Reading has been Finished" << std::endl ;
    gomfile.close() ; 
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    // Finding the centroid of each element
    #pragma omp parallel for schedule(static) shared(mesh,vertices_array) num_threads(8)
    for (size_t i = 0 ; i < mesh.size() ; ++i){

        mesh[i].centroidx = (vertices_array[mesh[i].v1].x + vertices_array[mesh[i].v2].x + vertices_array[mesh[i].v3].x) / 3 ; 
        mesh[i].centroidy = (vertices_array[mesh[i].v1].y + vertices_array[mesh[i].v2].y + vertices_array[mesh[i].v3].y) / 3 ;
    }
    
    // Bounding Box of each element 
    box.box_max_x = 1.0 ;   
    box.box_max_y = 1.0 ;
    box.box_min_x = -1.0 ;
    box.box_min_y = -1.0 ;
    box.centre_x = ( box.box_max_x + box.box_min_x ) * 0.5 ;
    box.centre_y = ( box.box_max_y + box.box_min_y ) * 0.5 ;
    
    //uint64_t N_level = 20 ; 
    uint64_t index = 0 ; 
    size_t meshsize = mesh.size() ; 

    // Generating Space Filling Curve 
    #pragma omp parallel for schedule(static) private(index,tempBox) shared(mesh,meshsize) num_threads(8)
    for (size_t i = 0 ; i < meshsize ; ++i) {
        
        tempBox.box_max_x = 1.0 ; 
        tempBox.box_max_y = 1.0 ; 
        tempBox.box_min_x = -1.0 ; 
        tempBox.box_min_y = -1.0 ;
        tempBox.centre_x = ( box.box_max_x + box.box_min_x ) * 0.5 ;
        tempBox.centre_y = ( box.box_max_y + box.box_min_y ) * 0.5 ;

        //mesh[i].sfc_index = Morton_Curve(mesh[i], box, 20) ;

        index = 0 ;
        for (size_t j = 0 ; j < 20 ; ++j) {
        
            index = index << 3;
        
            // Set the octant using Morton order curve
            if (mesh[i].centroidx > tempBox.centre_x) {
                index = index + 1 ;
            }
            if (mesh[i].centroidy > tempBox.centre_y) {
                index = index + 2 ;
            }
        
            // Update the bounding box
            if (mesh[i].centroidx > tempBox.centre_x) {
                tempBox.box_min_x = tempBox.centre_x ;
            }
            else tempBox.box_max_x = tempBox.centre_x ;
        
            if (mesh[i].centroidy > tempBox.centre_y) {
                tempBox.box_min_y = tempBox.centre_y ;
            }
            else tempBox.box_max_y = tempBox.centre_y ;
        
            // Update the bounding Box
            //updateBox(box) ;
            tempBox.centre_x = (tempBox.box_max_x + tempBox.box_min_x ) * 0.5 ;
            tempBox.centre_y = (tempBox.box_max_y + tempBox.box_min_y ) * 0.5 ;
        }
    
        mesh[i].sfc_index = index ; 

    }

    // Sorting of Elements 
    tbb::parallel_sort(mesh.begin(), mesh.end(), compareElement) ; // Quicksort of O(NlogN)
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout<< "Parallel SFC Generation took "<<std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()<< std::endl ; 
      
    // Renumbering the vertices according to space filling curves 
    for( auto it = mesh.begin() ; it!= mesh.end() ; ++it ){
             
            if(!vertices_array[it->v1].flag){
                vertices_array[it->v1].flag = true ; 
                vertices_array[it->v1].sfc_index = _globalVerticesNumber ; 
                ++_globalVerticesNumber ; 
            }
            
            if(!vertices_array[it->v2].flag){
                vertices_array[it->v2].flag = true ; 
                vertices_array[it->v2].sfc_index = _globalVerticesNumber ; 
                ++_globalVerticesNumber ; 
            }
                        
            if(!vertices_array[it->v3].flag){
                 
                vertices_array[it->v3].flag = true ; 
                vertices_array[it->v3].sfc_index = _globalVerticesNumber ; 
                ++_globalVerticesNumber ;
            }
    }
    
    

    int64_t tempindex = 0 ; 
    // Renumbering the vertices 
    #pragma omp parallel for  schedule(static) private(tempindex) shared(mesh,vertices_array) num_threads(8)
    for (size_t i = 0 ; i < mesh.size() ; ++i ){

        tempindex = vertices_array[mesh[i].v1].sfc_index ; 
        mesh[i].v1 = tempindex ; 

        tempindex = vertices_array[mesh[i].v2].sfc_index ; 
        mesh[i].v2 = tempindex ; 

        tempindex = vertices_array[mesh[i].v3].sfc_index ; 
        mesh[i].v3 = tempindex ; 

    }
    
    tbb::parallel_sort(vertices_array.begin(), vertices_array.end(), compareVertices) ; // QuickSort
 
    // // Rewriting the file 
    // std::ofstream outfile ; 
    // outfile.open("unit_circle_sfc.txt") ; 
    
    // outfile << _numVertices <<"\t vertices in the domain" << std::endl ; 
    // outfile << "Index of vertex | x cordinate | y cordiate "<< std::endl ; 
    
    // // Writing the vertices 
    // for (uint64_t i =0 ; i < _numVertices ; ++i){
            
    //     outfile <<"\t"<< i <<"\t"<<vertices_array[i].x<<"\t"<<vertices_array[i].y<< std::endl ; 
    // }
    
    
    // outfile<< _numElements <<" \t elements in the domain " << std::endl ; 
    // outfile<<"\tindex of vertex 0 | index of vertex 1 | index of vertex 2 "<< std::endl ;
    // // Writing the elements 
    // for (uint64_t i =0 ; i < mesh.size(); ++i ){
        
    //     outfile <<"\t"<<mesh[i].v1 <<"\t "<<mesh[i].v2<<"\t"<<mesh[i].v3 << std::endl ; 
        
    // }
    
    // outfile.close() ; 

}

bool compareVertices(const vertices &v1, const vertices &v2) {
    
    if (v1.sfc_index < v2.sfc_index) {
        return  true ;
    }
    else return false ;
}

bool compareElement(const element &e1, const element &e2) {
    
    if (e1.sfc_index <= e2.sfc_index) {
        return true ;
    }
    else
        return false;
}

void updateBox(BoundingBox &b) {
    
    b.centre_x = (b.box_max_x + b.box_min_x ) * 0.5 ;
    b.centre_y = (b.box_max_y + b.box_min_y ) * 0.5 ;

}

uint64_t Morton_Curve(element E,BoundingBox box,uint64_t N_level) {
    
    int64_t index = 0 ;
    //std::cout<<E.centroidx<<std::endl ;
    for (size_t i =0 ; i < N_level ; ++i) {
        
        index = index << 3;
        
        // Set the octant using Morton order curve
        if (E.centroidx > box.centre_x) {
            index = index + 1 ;
        }
        if (E.centroidy > box.centre_y) {
            index = index + 2 ;
        }
        
        // Update the bounding box
        if (E.centroidx > box.centre_x) {
            box.box_min_x = box.centre_x ;
        }
        else box.box_max_x = box.centre_x ;
        
        if (E.centroidy > box.centre_y) {
            box.box_min_y = box.centre_y ;
        }
        else box.box_max_y = box.centre_y ;
        
        // Update the bounding Box
        updateBox(box) ;
    }
    return  index ;
}












