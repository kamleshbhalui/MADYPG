#ifndef FBXNODE_H
#define FBXNODE_H

#include "fbxproperty.h"
#include <algorithm>

namespace fbx {

class FBXNode
{
public:
    FBXNode();
    FBXNode(std::string name);

    std::uint32_t read(std::ifstream &input, uint32_t start_offset);
    std::uint32_t write(std::ofstream &output, uint32_t start_offset);
    void print(std::string prefix="");
    bool isNull();

    void addProperty(int16_t);
    void addProperty(bool);
    void addProperty(int32_t);
    void addProperty(float);
    void addProperty(double);
    void addProperty(int64_t);
    void addProperty(const std::vector<bool>);
    void addProperty(const std::vector<int32_t>);
    void addProperty(const std::vector<float>);
    void addProperty(const std::vector<double>);
    void addProperty(const std::vector<int64_t>);
    void addProperty(const std::vector<uint8_t>, uint8_t type);
    void addProperty(const std::string);
    void addProperty(const char*);
    void addProperty(FBXProperty);

    void addPropertyNode(const std::string name, int16_t);
    void addPropertyNode(const std::string name, bool);
    void addPropertyNode(const std::string name, int32_t);
    void addPropertyNode(const std::string name, float);
    void addPropertyNode(const std::string name, double);
    void addPropertyNode(const std::string name, int64_t);
    void addPropertyNode(const std::string name, const std::vector<bool>);
    void addPropertyNode(const std::string name, const std::vector<int32_t>);
    void addPropertyNode(const std::string name, const std::vector<float>);
    void addPropertyNode(const std::string name, const std::vector<double>);
    void addPropertyNode(const std::string name, const std::vector<int64_t>);
    void addPropertyNode(const std::string name, const std::vector<uint8_t>, uint8_t type);
    void addPropertyNode(const std::string name, const std::string);
    void addPropertyNode(const std::string name, const char*);

    void addChild(FBXNode child);
    uint32_t getBytes();


    fbx::FBXNode* getChildByName(const std::string name) {
        for (fbx::FBXNode& child : children) {
            if (child.getName() == name) {
                return &child;
            }
        }
        return nullptr;
    }
    std::vector<fbx::FBXNode*> getChildrenByName(const std::string name) {
        std::vector<fbx::FBXNode*> result;
        for (fbx::FBXNode& child : children) {
            if (child.getName() == name) {
                result.push_back(&child);
            }
        }
        return result;
    }

    void overwriteProperty(FBXProperty && prop) {
        properties.clear();
        properties.push_back(prop);
    }

    void removeChild(const std::string& name) {
        children.erase(
            std::remove_if(children.begin(), children.end(),
                [&](const FBXNode & o) { return o.name == name; }),
            children.end());
    }


    const std::vector<FBXNode> getChildren();
    const std::string getName();
private:
    std::vector<FBXNode> children;
    std::vector<FBXProperty> properties;
    std::string name;
};

} // namespace fbx

#endif // FBXNODE_H
