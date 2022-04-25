#include "..\include\cell.h"

/** Implementation **/
Node Cell::centroid(const std::vector<Node> &nodes)
{
    Node c;
    Node v1 = nodes[nodeIndices[0]],
         v2 = nodes[nodeIndices[1]],
         v3 = nodes[nodeIndices[2]];
    c.x = (v1.x + v2.x + v3.x) / 3.;
    c.y = (v1.y + v2.y + v3.y) / 3.;
    return c;
}
Node Cell::centroid(const Node &v1, const Node &v2, const Node &v3)
{
    Node c;
    c.x = (v1.x + v2.x + v3.x) / 3.;
    c.y = (v1.y + v2.y + v3.y) / 3.;
    return c;
}

double Cell::volume(const std::vector<Node> &nodes)
{
    double vol_;
    Node v1 = nodes[nodeIndices[0]],
         v2 = nodes[nodeIndices[1]],
         v3 = nodes[nodeIndices[2]];
    vol_ = (v1.x - v2.x) * (v1.y + v2.y);
    vol_ += (v2.x - v3.x) * (v2.y + v3.y);
    vol_ += (v3.x - v1.x) * (v3.y + v1.y);
    return 0.5 * vol_;
}
double Cell::volume(const Node &v1, const Node &v2, const Node &v3)
{
    double vol_;
    vol_ = (v1.x - v2.x) * (v1.y + v2.y);
    vol_ += (v2.x - v3.x) * (v2.y + v3.y);
    vol_ += (v3.x - v1.x) * (v3.y + v1.y);
    return 0.5 * vol_;
}
double Cell::volume_quad(const Node &v1, const Node &v2, const Node &v3, const Node &v4)
{
    return volume(v1, v2, v3) + volume(v1, v3, v4);
}

std::ostream &operator<<(std::ostream &os, const Cell &cell)
{
    os << cell.index << '\t'; // cell.VTK_type << "\t";
    for (const auto &nIndx : cell.nodeIndices)
        os << nIndx << "\t";
    os << cell.boundaryFlag;
    return os;
}