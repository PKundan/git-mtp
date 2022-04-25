#include "..\include\readMesh.h"


/********************************************************************************
/ readCell Function : reads the line containing cell info(su2 format)           /
/ Inputs :- line : string , line containing cell info                           /
/           cells : vector of cells                                             /
/ Result :- cell read from the line will be appended in cells vector            /
/*******************************************************************************/

void readCell(const std::string &line,
              std::vector<Cell> &cells)
{
    Cell _cell;                 // temporary cell object
    std::vector<double> vals;   // vector of values to be read froom line
    std::stringstream ss(line); // creating stringstream of line

    while (ss.good()) // checks if it's not the end of stringstream
    {
        std::string substr;
        getline(ss, substr, ' '); // creating substrings seperated by space
        std::stringstream val(substr);
        double dummy_val;
        val >> dummy_val;
        // std::cout<< dummy_val << "\t";
        vals.push_back(dummy_val); // putting values read in vals vector
    }
    // _cell.VTK_type = vals.front(); // first value is VTK_type
    for (size_t j = 1; j < vals.size() - 1; ++j)
    {
        _cell.nodeIndices.push_back(vals[j]); // middle values are nodeIndices
    }
    _cell.index = vals.back(); // last value is cell-index
    cells.push_back(_cell);    // appending cell in cells vector
}

/********************************************************************************
/ readNode Function : reads the line containing node info(su2 format)           /
/ Inputs :- line : string , line containing node info                           /
/           nodes : vector of nodes                                             /
/ Result :- node read from the line will be appended in nodes vector            /
/*******************************************************************************/
void readNode(const std::string &line,
              std::vector<Node> &nodes)
{
    std::vector<double> vals;   // vector of values to be read froom line
    std::stringstream ss(line); // creating stringstream of line

    while (ss.good())
    {
        std::string substr;
        getline(ss, substr, ' '); // creating substrings seperated by space
        std::stringstream val(substr);
        double dummy_val;
        val >> dummy_val;
        // std::cout<< dummy_val << "\t";
        vals.push_back(dummy_val); // appending values in vals vector
    }
    nodes.push_back(Node(vals[0], vals[1])); // first value : x-cooordinate, second: y_coordinate
    //.........................................// appending node in  nodes vector
}

/********************************************************************************
/ readEdge Function : reads the line containing boundary edge info(su2 format)  /
/ Inputs :- line : string , line containing boundary edge info                  /
/           maekerTag : tag of the boundary, ex. inlet, outlet, wall            /
/           boundaryEdges : vector of edges at the boundary                     /
/ Result :- edges read from the line will be appended in boundaryEdges vector   /
/*******************************************************************************/
void readEdge(std::string &line,
              const std::string &markerTag,
              std::vector<Edge> &boundaryEdges)
{
    int flag;
    if (markerTag == "inlet")       // for setting boundaryFlag
        flag = 1;                   // inelt : 1
    else if (markerTag == "outlet") // outlet : 2
        flag = 2;                   // wall : 3
    else if (markerTag == "wall")
        flag = 3;
    else
        flag = 3;

    std::vector<double> vals;   // vector of values to be read froom line
    std::stringstream ss(line); // creating stringstream of line

    while (ss.good())
    {
        std::string substr;
        getline(ss, substr, ' '); // creating substrings seperated by space
        std::stringstream val(substr);
        double dummy_val;
        val >> dummy_val;
        vals.push_back(dummy_val); // appending values in vals vector
    }
    Edge edg;                                    // creating temporary Edge object
    edg.index = boundaryEdges.size();            // initially size of boundaryEdges will be zero. //It will increase as edges are appended
    for (size_t k = 1; k < vals.size() - 1; ++k) // It will increase as edges will get appended in it.
        edg.nodeIndices.push_back(vals[k]);      // middle values are nodeIndices
    edg.boundaryFlag = flag;                     // setting boundaryFlag
    boundaryEdges.push_back(edg);                // appending edge in boundaryEdges
}
/*******************************************************************************
/ readMesh Function : Modifies input vectors nodes, cells and boundaryedges     /
/                       by reading the meshfile( *.su2 format)                  /
/ Inputs :- fileName : string , mesh file name (*.su2 format)                   /
/           nodes : vector of <Node> objects, size = 0; (initially empty)       /
/           cells : vector of <Cell> objects, size = 0; (initially empty)       /
/           boundaryEdges : vector of <Edge> objects, size=0, (initailly empty) /
/ Result :-  modification in nodes, cells and boundaryEdges vectors             /
/*******************************************************************************/
void readMesh(std::string &fileName,
              std::vector<Node> &nodes,
              std::vector<Cell> &cells,
              std::vector<Edge> &boundaryEdges)
{
    std::ifstream fin(fileName, std::ifstream::binary); // Open file to read
    if (fin)
    {
        while (!fin.eof()) // if not end of file
        {
            std::string line;
            getline(fin, line);                 // reading line
            if (line.find('=') < line.length()) // find if line contains '='
            {
                line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end()); // remove whitespaces in line
                auto delimiterPos = line.find("=");                                        // find position of '='
                auto name = line.substr(0, delimiterPos);                                  // name = substring before '='
                auto value = line.substr(delimiterPos + 1);                                // value  = substring after '='
                int int_val;
                std::stringstream sVal(value);
                sVal >> int_val; // converting value from 'string' to 'int'
                // std::cout << name << "=" << value << std::endl;
                if (name == "NELEM") // NELEM = np. of elements or cells
                {
                    for (int i = 0; i != int_val; ++i) // int_val = no. of cells/elements, iterating
                    {                                  // reading cells
                        getline(fin, line);            // reading line
                        readCell(line, cells);         // updating cell info and pushing cell into cells vector
                    }
                }
                else if (name == "NPOIN")  // NPOIN : no. of nodes
                { // Reading nodes
                    for (int i = 0; i != int_val; ++i) // int_val = no. of nodes/points, iterating
                    {
                        getline(fin, line);
                        readNode(line, nodes);  // reading node info and pushing node into nodes vector
                    }
                }
                else if (name == "NMARK") // NMARK : no. of boundary markers, ex. inlet, outlet, wall
                {
                    for (int i = 0; i != int_val; ++i) // int_val = NMARK
                    {
                        getline(fin, line);
                        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end()); //remove whitespaces
                        delimiterPos = line.find("=");  // find position of '=' delimiter
                        auto markerTag = line.substr(delimiterPos + 1); // markerTag = subtring after '=', ex. inlet, wall
                        // std::cout << markerTag << std::endl;
                        getline(fin, line);
                        line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end()); //remove whitespaces
                        delimiterPos = line.find("=");  // find position of '=' delimiter
                        value = line.substr(delimiterPos + 1); // value = subtring after '=', no. of marker elements, ex. no. of inlet boundary edges 
                        std::stringstream ss(value);
                        int nMarkElems;
                        ss >> nMarkElems;   // no. of marker elements (int)
                        for (int k = 0; k != nMarkElems; ++k) // iterating to read all edges
                        {
                            getline(fin, line);
                            readEdge(line, markerTag, boundaryEdges); // read edge info and append in boundaryEdges
                        }
                    }
                }
            }
        }
    }
    else
    {
        std::cout << "File does not exist!" << std::endl;
    }
}