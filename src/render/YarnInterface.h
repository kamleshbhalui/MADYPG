#ifndef __YARNINTERFACE__H__
#define __YARNINTERFACE__H__

#include "EigenDefinitions.h"

class YarnInterface {
    public:
        YarnInterface() {}
        virtual ~YarnInterface() {}

        virtual void step() {}
        virtual const std::vector<uint32_t>& getIndices() = 0;
        // virtual const MatrixXXs& getVertexData() = 0;
};

#endif // __YARNINTERFACE__H__