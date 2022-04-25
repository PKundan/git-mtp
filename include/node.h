#ifndef _NODE_H
#define _NODE_H

#include<iostream>

/*****************************************************************/
// ----------------------Node class-------------------------------/
/*****************************************************************/
/* Attributes:
    x : x co-ordinate
    y : y co-ordinate

/* Operators :
        =    :  copy assignment operator
        <<   :  to print node attributes on console
/******************************************************************/

class Node{
public:
    double x, y;
    int flag = 0;  // vertex = 0; cell_centroid = 1; edge_center = 2;
    Node(const double& x_ = 0, const double& y_ = 0) : x(x_), y(y_) {}
    // Copy constructor
    Node(const Node& p1); //{ x = p1.x; y = p1.y; }
    Node operator=(const Node& v); 
    Node& operator+=(const Node& other);
    friend std::ostream& operator<<(std::ostream& os, const Node& node);
};

#endif



