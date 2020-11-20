#ifndef FBXDOCUMENT_H
#define FBXDOCUMENT_H

#include "fbxnode.h"

namespace fbx {

class FBXDocument
{
public:
    FBXDocument();
    void read(std::ifstream &input);
    void read(std::string fname);
    void write(std::string fname);
    void write(std::ofstream &output);

    void createBasicStructure();

    std::vector<FBXNode> nodes;


    fbx::FBXNode* getNodeByName(const std::string name) {
        for (fbx::FBXNode& child : nodes) {
            if (child.getName() == name) {
                return &child;
            }
        }
        return nullptr;
    }

    std::uint32_t getVersion();
    void print();

private:
    std::uint32_t version;
};

} // namespace fbx

#endif // FBXDOCUMENT_H
