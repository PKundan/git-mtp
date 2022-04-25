#include "..\include\node.h"

// Node::Node(const double& x_ = 0, const double& y_ = 0) : x(x_), y(y_) {}
Node::Node(const Node &p1)
{
    x = p1.x;
    y = p1.y;
    flag = p1.flag;
}

Node Node::operator=(const Node &v)
{
    this->x = v.x;
    this->y = v.y;
    this->flag = v.flag;
    return *this;
}

Node& Node::operator+=(const Node& other){
    x += other.x;
    y += other.y;
    return *this;
}

Node operator+(Node lhs, const Node& rhs){
    lhs += rhs;
    return lhs;
}

Node operator*(Node n, const double& k){
    n.x *= k;
    n.y *= k;
    return n;
}

Node operator*(const double& k, Node n){
    return n * k;
}


std::ostream& operator<<(std::ostream& os, const Node& node)
{
    os << node.x << '\t' << node.y << '\t' << node.flag;
    return os;
}
