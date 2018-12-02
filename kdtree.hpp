/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>
#include <stack>

using namespace std;

template <int Dim>
double KDTree<Dim>::distance(const Point<Dim>& first, const Point<Dim>& second) const
{
    double d = 0;
    for (int i = 0; i < Dim; i++) {
        d += (first[i] - second[i]) * (first[i] - second[i]);
    }
    return d;
}

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
    if (first[curDim] == second[curDim])
        return first < second;
    else
        return first[curDim] < second[curDim];
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
    double d1 = distance(target, currentBest);
    double d2 = distance(target, potential);
    if (d2 == d1)
        return potential < currentBest;
    else
        return d2 < d1;
}

template <int Dim>
int KDTree<Dim>::quickPartition(vector<Point<Dim>>& points, int dim, int left, int right, int pivotIndex)
{
    Point<Dim> pivotValue = points[pivotIndex];

  
    swap(points[pivotIndex], points[right]);

    int storeIndex = left;
    for (int i = left; i < right; i++) {
        if (smallerDimVal(points[i], pivotValue, dim)) {
            swap(points[storeIndex], points[i]);
            
            storeIndex++;
        }
    }
    swap(points[storeIndex], points[right]);

    return storeIndex;
}

template <int Dim>
int KDTree<Dim>::quickSelect(vector<Point<Dim>>& points, int dim, int left, int right, int k)
{
    while(1){
        if (left == right)
            return left;
        int pivotIndex = (left + right) / 2;
        pivotIndex = quickPartition(points, dim, left, right, pivotIndex);
        if (k == pivotIndex)
            return k;
        else if (k < pivotIndex)
            right = pivotIndex - 1;
        else
            left = pivotIndex + 1;
    }
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode *KDTree<Dim>::build(vector<Point<Dim>> &points, int dim, int start, int end) {
    int pivot = (start + end) / 2;
    if(end - start < 0){
        return NULL;
    }


    pivot = quickSelect(points, dim, start, end, pivot);

    KDTreeNode *node = new KDTreeNode(points[pivot]);
        node->left = build(points, (dim + 1) % Dim, start, pivot - 1);
   
        node->right = build(points, (dim + 1) % Dim, pivot + 1, end);
    return node;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
   
    if (newPoints.size() == 0) {
        root = nullptr;
        size = 0;
        return;
    }

    vector<Point<Dim>> points;
    for (int i = 0; i < (int)newPoints.size(); i++) {
        points.push_back(newPoints[i]);
    }

    root = build(points, 0, 0, points.size()-1);
    size = points.size();
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  root = copy(other.root);
  

// Copy this node and its children
 
}


template <int Dim>
typename KDTree<Dim>::KDTreeNode * KDTree<Dim>::copy(const KDTreeNode * subRoot)
{
    if (subRoot == NULL)
        return NULL;

    // Copy this node and its children
    KDTreeNode* newNode = new KDTreeNode(subRoot->point);
    newNode->left = copy(subRoot->left);
    newNode->right = copy(subRoot->right);
    return newNode;
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   * 
   */
  if (this != &rhs) {
    clear(root);
    root = copy(rhs.root);
}
    size = rhs.size();

  return *this;
}

// template <int Dim>
// KDTree<Dim>::clear()

template <int Dim>
void KDTree<Dim>::clear(KDTreeNode* subRoot)
{
    if (subRoot == NULL)
        return;
    clear(subRoot->left);
    clear(subRoot->right);
    delete subRoot;
}


template <int Dim>
KDTree<Dim>::~KDTree() {
    /**
     * @todo Implement this function!
     */

    std::stack<KDTreeNode*> nodes;
    if (root != nullptr)
        nodes.push(root);
    while (!nodes.empty()) {
        KDTreeNode *node = nodes.top();
        nodes.pop();
        if (node->left)
            nodes.push(node->left);
        if (node->right)
            nodes.push(node->right);
        delete node;
    }

    root = nullptr;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const KDTreeNode *node, const Point<Dim>& query, const int dim) const
{
    Point<Dim> it = node->point;
    Point<Dim> x;
    double radiusSquared = distance(query, it);
    int dir = 0;

    if (node->left != nullptr && smallerDimVal(query, node->point, dim)) {
        dir = 0;
        x = findNearestNeighbor(node->left, query, (dim + 1) % Dim);
    } else if (node->right != nullptr && smallerDimVal(node->point, query, dim)) {
        dir = 1;
        x = findNearestNeighbor(node->right, query, (dim + 1) % Dim);
    } else {
        x = it;
    }

    if (shouldReplace(query, node->point, x)) {
        it = x;
    } else {
        it = node->point;
    }
    radiusSquared = distance(query, it);

    double planeDistanceSquared = abs(query[dim] - node->point[dim]);
    planeDistanceSquared = planeDistanceSquared * planeDistanceSquared;

    if (planeDistanceSquared <= radiusSquared) {
        if (dir == 0 && node->right != nullptr) {
            x = findNearestNeighbor(node->right, query, (dim + 1) % Dim);
            if (shouldReplace(query, it, x)) {
                it = x;
            }
        }
        if (dir != 0 && node->left != nullptr) {
            x = findNearestNeighbor(node->left, query, (dim + 1) % Dim);
            if (shouldReplace(query, it, x)) {
                it = x;
            }
        }
    }

    return it;
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    return findNearestNeighbor(root, query, 0);
}


// /**
//  * @file kdtree.cpp
//  * Implementation of KDTree class.
//  */

// #include <utility>
// #include <algorithm>

// using namespace std;

// template <int Dim>
// bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
//                                 const Point<Dim>& second, int curDim) const
// {
//     /**
//      * @todo Implement this function!
//      */
//     if(first[curDim] < second[curDim]){
//         return true;
//     }else if(first[curDim] == second[curDim]){
//         return first < second;
//     }

//     return false;
// }

// template <int Dim>
// bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
//                                 const Point<Dim>& currentBest,
//                                 const Point<Dim>& potential) const
// {
//     /**
//      * @todo Implement this function!
//      */
//     double total = 0;
//     double newTotal = 0;
//     for(unsigned int i = 0; i < Dim; i++){
//         total = total + (target[i] - currentBest[i]) * (target[i] - currentBest[i]);
//     }
    
//     for(unsigned int i = 0; i < Dim; i++){
//         newTotal = newTotal + (target[i] - potential[i]) * (target[i] - potential[i]);
//     }

//     if(newTotal < total){
//         return true;
//     }

//      return false;
// }

// template <int Dim>
// KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
// {
//     /**
//      * @todo Implement this function!
//      * 
//      */
//     if(newPoints.size() == 1){
//         root = new KDTreeNode(newPoints[0]);

//     }
//     int size = newPoints.size();

//     createTree(newPoints, 0, size-1, 0, root);
    
// }


// template <int Dim>
// int KDTree<Dim>::partition(const vector<Point<Dim>>& newPoints, int pivotPoint, int left, int right, int dim){
//     swap(newPoints[pivotPoint], newPoints[right]);
//     int storeIndex = left; 
//     for(unsigned int i = left; i < right - 1; i++){
//         if(smallerDimVal(newPoints[i], newPoints[pivotPoint], dim)){
//             swap(newPoints[storeIndex], newPoints[i]);
//             storeIndex++;
//         } 
//     }
//     swap(newPoints[storeIndex], newPoints[right]);
//     return storeIndex;
// }

// template <int Dim>
// Point<Dim> KDTree<Dim>::select(const vector<Point<Dim>>& newPoints, int left, int right, int dim, int k){
//     if(left == right){
//         return newPoints[left];
//     }
    
//     int pivotPoint = (left+right)/2;
//     pivotPoint = partition(newPoints, pivotPoint, left, right, dim);
//     if(k == pivotPoint){
//         return newPoints[k];
//     }else if(k < pivotPoint){
//         return select(newPoints, left, pivotPoint - 1, dim, k);
//     }else{
//         return select(newPoints, pivotPoint + 1, left, dim, k);
//     }
    
// }
// template <int Dim>
// void KDTree<Dim>::createTree(const vector<Point<Dim>>& newPoints,int left, int right, int dim, KDTreeNode *& cur){
//     if(right>left){
//         return;
//     }
    
//     cur = new KDTreeNode(select(newPoints, left, right, dim, (left+right)/2));
//   createTree(newPoints, left, right, (dim+1)%Dim, cur->left);
//   createTree(newPoints, left, right, (dim=1)%Dim, cur->right);
    
// }


// template <int Dim>
// KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
//   /**
//    * @todo Implement this function!
//    */
// }

// template <int Dim>
// const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
//   /**
//    * @todo Implement this function!
//    */

//   return *this;
// }

// template <int Dim>
// KDTree<Dim>::~KDTree() {
//   /**
//    * @todo Implement this function!
//    */
// }

// template <int Dim>
// Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
// {
//     /**
//      * @todo Implement this function!
//      */

//     return Point<Dim>();
// }

