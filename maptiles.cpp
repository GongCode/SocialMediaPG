/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
#include "kdtree.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    MosaicCanvas * newCanvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());

    vector<Point<3>> colorVector;
    for(unsigned int c = 0; c < theTiles.size(); c++){
        LUVAPixel pixel = theTiles[c].getAverageColor();
        colorVector.push_back(convertToXYZ(pixel));
    }
    KDTree<3> tree(colorVector);

    

    for(int i = 0; i < newCanvas->getRows(); i++){
        for(int j = 0; j < newCanvas->getColumns(); j++){
            LUVAPixel canvasPixel = theSource.getRegionColor(i,j);
            Point<3> currentPixel = convertToXYZ(canvasPixel);
            Point<3> nearestPixel = tree.findNearestNeighbor(currentPixel); 
            unsigned int z;
            for(z = 0; z < theTiles.size(); z++){
                if(nearestPixel == colorVector[z]){
                    break;
                }
            }
            newCanvas->setTile(i,j, &theTiles[z]);
        }

    }
    


    return newCanvas;
}

