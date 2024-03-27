#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <complex>
#include <vector>
#include <cmath>
#include <Eigen/Eigenvalues>

// #define INTERACTIONS
// #define PROXIMITY

static const std::complex< double > iDelta(0, 0.01);
static const double omegaMin = -3.5;
static const double omegaMax =  12.5;
static const size_t omegaPoints = 10001;

static const double PI = 4.0 * atan(1.0);

class Node
{
public:
    double J;
    Node(double coupling) {
        J = coupling;
        m_energy = 2 * J;
    }
    void setValue(std::complex< int > value) { m_value = value; };
    void setMultiplicity(short multiplicity) { m_multiplicity = multiplicity; };
    void setEnergy(double energy) { m_energy = energy; };
    void setParent(Node *parent) { m_parent = parent; };
    void pushChild(size_t id, std::complex< int > value, short multiplicity) {
        m_child[id] = new Node(J);
        m_child[id]->setParent(this);
        m_child[id]->setValue(value);
        m_child[id]->setMultiplicity(multiplicity);
        m_child[id]->setEnergy(this->getEnergy() + 2 * J);
    };
    std::complex< int > getValue() { return m_value; };
    short getMultiplicity() { return m_multiplicity; };
    double getEnergy() { return m_energy; };
    Node *getParent() { return m_parent; };
    Node *getChild(size_t id) { return m_child[id]; };
private:
    std::complex< int > m_value = 0;
    short m_multiplicity = 1;
    double m_energy;
    Node *m_parent = NULL;
    Node *m_child[4] = {};
};

template <typename T>
void printToFile(T output, std::string fileName)
{
    std::ofstream file;
    file.open(fileName);
    if(file.good()) {
        file << output;
        file.close();
    }
};

using namespace Eigen;

bool init(int, char const *[], short &, double &);
void build(Node *, short);
void span(Node *, short);
void print(Node *, short);
void scanPrint(Node *, std::vector< std::complex< int > >, short, short);
void prepare(Node *, short, double J, bool, bool);
size_t systemInfo(Node *, short, bool);
size_t getSystemInfo(Node *, short, std::vector< size_t > &);
void searchGraph(Node *, size_t &, std::vector< size_t > &, short, short);
void calculateEnergy(Node *, short, double);
void calculateEnergyRecursively(Node *, short, short, double);
std::vector< std::string > readPaths(std::string);
std::string getFileName(std::string, double);
Node *getNode(Node *, std::string);
void calculate(Node *, std::string, short, double);
std::complex< double > getSelfEnergyValue(double, Node *);
VectorXcd calculateSelfEnergy(Node *);

int main(int argc, char const *argv[])
{
    short maxLength = 0;
    double J = 0.4;

    if(!init(argc, argv, maxLength, J)) return 1;

    Node *graph = new Node(J);
    build(graph, maxLength);
    prepare(graph, maxLength, J, false, false); // (shouldShowWalks, shouldShowStatesCount)
    auto paths = readPaths("input.txt");
    for (auto path : paths) {
        calculate(graph, path, maxLength, J);
    }
    return 0;
}

bool init(int argc, char const *argv[], short &maxLength, double &J)
{
    if(argc > 1) {
        std::stringstream s(argv[1]);
        s >> maxLength;
    }
    if(argc > 2) {
        std::stringstream s(argv[2]);
        s >> J;
    }
    if(maxLength < 0 || maxLength > 20) {
        std::cout << "Requested Number of Magnons is Negative or Too High" << std::endl;
        return false;
    }
    return true;
}

void build(Node *graph, short maxPathLength)
{
    std::cout << "Building Walks Graph - IN PROGRESS...";
    if(maxPathLength > 0) {
        graph->pushChild(0, std::complex< int >(1, 0), 4);
        span(graph->getChild(0), maxPathLength - 1);
    }
    std::cout << "\n" << "Building Walks Graph - FINISHED      " << std::endl;
}

void span(Node *node, short remainingSpans)
{
    if(remainingSpans > 0) {
        if(node->getMultiplicity() == 4) {
            node->pushChild(0, node->getValue() + std::complex< int >(1, 0), 4);
            span(node->getChild(0), remainingSpans - 1);
            node->pushChild(1, node->getValue() + std::complex< int >(0, 1), 8);
            span(node->getChild(1), remainingSpans - 1);
        } else {
            using namespace std::complex_literals;
            std::vector< std::complex< int > > move = {1, 1i, -1, -1i};
            for(auto it = 0; it < 4; it++) {
                bool isNewNode = true;
                std::complex< int > newValue = node->getValue() + move[it];
                Node *parent = node->getParent();
                while(parent) {
                    if(parent->getValue() == newValue) {
                        isNewNode = false;
                        break;
                    }
                    parent = parent->getParent();
                }
                if(isNewNode) {
                    node->pushChild(it, newValue, 8);
                    span(node->getChild(it), remainingSpans - 1);
                }
            }
        }
    }
}

void prepare(Node *graph, short maxLength, double J, bool shouldShowWalks, bool shouldShowStatesCount)
{
    size_t basisSize = systemInfo(graph, maxLength, shouldShowStatesCount);
    calculateEnergy(graph, maxLength, J);
    if(shouldShowWalks) print(graph, maxLength);
}

void print(Node *graph, short maxDepth)
{
    std::cout << "Representative Self-avoiding Walks up to " << maxDepth << " Magnons -- START" << std::endl;
    std::vector< std::complex< int > > path = {graph->getValue()};
    std::cout << " (" << 0 << ": " << path[0] << ")" << " | Energy : " << graph->getEnergy() << std::endl;
    if(maxDepth > 0)
        scanPrint(graph, path, 0, maxDepth);
    std::cout << "Representative Self-avoiding Walks up to " << maxDepth << " Magnons -- END" << std::endl;
}

void scanPrint(Node *parent, std::vector< std::complex< int > > path, short depth, short maxDepth)
{
    if(depth < maxDepth) {
        for(auto it = 0; it < 4; it++) {
            Node *node = parent->getChild(it);
            if(node) {
                path.push_back(node->getValue());
                std::cout << " (" << 0 << ": " << path[0] << ")";
                for(auto it = 1; it < path.size(); it++)
                    std::cout << " -> " << "(" << it << ": " << path[it] << ")";
                std::cout << " | Energy : " << node->getEnergy() << std::endl;
                scanPrint(node, path, depth + 1, maxDepth);
                path.erase(path.begin() + path.size() - 1);
            }
        }
    }
}

size_t systemInfo(Node *graph, short maxLength, bool shouldShowStatesCount)
{
    std::vector< size_t > selfAvoidingWalksCount;
    size_t basisSize = getSystemInfo(graph, maxLength, selfAvoidingWalksCount);
    std::cout << "System Info -- START" << std::endl;
    if(shouldShowStatesCount) {
        std::cout << " (Steps -> States):" << std::endl;
        for(auto it = 0; it < selfAvoidingWalksCount.size(); it++)
            std::cout << " (" << it << " -> " << selfAvoidingWalksCount[it] << ")" << std::endl;
    }
    std::cout << " Basis Size: " << basisSize << std::endl;
    std::cout << "System Info -- END" << std::endl;
    return basisSize;
}

size_t getSystemInfo(Node *graph, short maxPathLength, std::vector< size_t > &selfAvoidingWalksCount)
{
    selfAvoidingWalksCount.push_back(1);
    size_t basisSize = 1;
    if(maxPathLength > 0)
        searchGraph(graph, basisSize, selfAvoidingWalksCount, 1, maxPathLength);
    return basisSize;
}

void searchGraph(Node *parent, size_t &basisSize, std::vector< size_t > &selfAvoidingWalksCount, short pathLength, short maxPathLength)
{
    if(pathLength <= maxPathLength) {
        for(auto it = 0; it < 4; it++) {
            Node *node = parent->getChild(it);
            if(node) {
                if(node->getMultiplicity() == 4) selfAvoidingWalksCount.push_back(4);
                else selfAvoidingWalksCount[pathLength] += 8;
                searchGraph(node, ++basisSize, selfAvoidingWalksCount, pathLength + 1, maxPathLength);
            }
        }
    }
}

void calculateEnergy(Node *graph, short maxPathLength, double J)
{
    std::cout << "Calculating Energies - IN PROGRESS...";
    if(maxPathLength > 0) {
        calculateEnergyRecursively(graph, 1, maxPathLength, J);
    }
    std::cout << "\n" << "Calculating Energies - FINISHED      " << std::endl;
}

#ifdef INTERACTIONS
void calculateEnergyRecursively(Node *parent, short pathLength, short maxPathLength, double J)
{
    short holeNeighboursCount = 0;
    if(pathLength <= maxPathLength) {
        holeNeighboursCount = 4;
        for(auto it = 0; it < 4; it++) {
            Node *node = parent->getChild(it);
            if(node)
                holeNeighboursCount -= node->getMultiplicity() / parent->getMultiplicity();
        }
        for(auto it = 0; it < 4; it++) {
            Node *node = parent->getChild(it);
            if(node) {
                node->setEnergy(parent->getEnergy() + J * (2.0 - static_cast< double >(holeNeighboursCount)));
                calculateEnergyRecursively(node, pathLength + 1, maxPathLength, J);
            }
        }
    }
#else
void calculateEnergyRecursively(Node *parent, short pathLength, short maxPathLength, double J)
{
#ifdef PROXIMITY
    short holeNeighboursCount = 0;
#endif
    if(pathLength <= maxPathLength) {
#ifdef PROXIMITY
        holeNeighboursCount = 4;
#endif
        for(auto it = 0; it < 4; it++) {
            Node *node = parent->getChild(it);
            if(node) {
#ifdef PROXIMITY
                holeNeighboursCount -= node->getMultiplicity() / parent->getMultiplicity();
#endif
                calculateEnergyRecursively(node, pathLength + 1, maxPathLength, J);
            }
        }
    }
#endif
#ifdef PROXIMITY
    else {
        using namespace std::complex_literals;
        std::vector< std::complex< int > > move = {1, 1i, -1, -1i};
        for(auto it = 0; it < 4; it++) {
            std::complex< int > value = parent->getValue() + move[it];
            Node *node = parent->getParent();
            while(node) {
                if(node->getValue() == value)
                    holeNeighboursCount++;
                node = node->getParent();
            }
        }
    }
    parent->setEnergy(parent->getEnergy() - 0.5 * J * static_cast< double >(holeNeighboursCount));
#endif
}

std::vector< std::string > readPaths(std::string fileName)
{
    std::ifstream file;
    std::vector< std::string > paths;
    file.open(fileName);
    if(file.good()) {
        std::string line;
        while(getline(file, line)) {
            paths.push_back(line);
        }
        file.close();
    }
    return paths;
}

std::string getFileName(std::string path, double J)
{
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4) << J;

    std::string fileName = "output/";
    fileName += "/J=" + ss.str() + "/SE_" + path + ".txt";

    return fileName;
}

Node *getNode(Node *graph, std::string path)
{
    std::vector< char > D = {'E', 'N', 'W', 'S'};
    for (size_t it = 0; it < path.size(); it++) {
        for (size_t id = 0; id < 4; id++) {
            if (path[it] == D[id]) {
                graph = graph->getChild(id);
                break;
            }
        }
    }
    return graph;
}

void calculate(Node *graph, std::string path, short maxLength, double J)
{
    std::cout << "Calculating Self-Energy " << path << " - IN PROGRESS...";
    VectorXcd result = calculateSelfEnergy(getNode(graph, path));
    IOFormat precise(FullPrecision, DontAlignCols, ", ", "\n");
    printToFile(result.format(precise), getFileName(path, J));
    std::cout << "\n" << "Calculating Self-Energy " << path << " - FINISHED      " << std::endl;
}

std::complex< double > getSelfEnergyValue(double omega, Node *node)
{
    std::complex< double > result = 0.0;
    if(node) {
        for(auto it = 0; it < 4; it++) {
            Node *child = node->getChild(it);
            if(child)
                result += getSelfEnergyValue(omega, child);
        }
        result = static_cast< std::complex< double >>(node->getMultiplicity() / node->getParent()->getMultiplicity()) / (omega + iDelta - node->getEnergy() - result);
    }
    return result;
}

VectorXcd calculateSelfEnergy(Node *graph)
{
    VectorXcd result = VectorXcd::Zero(omegaPoints);
    #pragma omp parallel for
    for(auto it = 0; it < omegaPoints; it++) {
        double omega = omegaMin + static_cast< double >(it) * (omegaMax - omegaMin) / static_cast< double >(omegaPoints - 1);
        result(it) = getSelfEnergyValue(omega, graph);
    }
    return result;
}
