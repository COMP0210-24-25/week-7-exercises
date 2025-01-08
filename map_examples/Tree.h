#include <stdexcept>
#include <memory>
#include <ostream>

template<typename keyType, typename valType>
class BinaryTree
{
    public:
    BinaryTree(keyType k, valType v) : key(k), value(v)
    {
        parent = nullptr;
    }

    void insert(const keyType k, const valType v)
    {
        if(k == key)
        {
            value = v;
        }
        else if(k < key)
        {
            insertChild(k, v, left);
        }
        else
        {
            insertChild(k, v, right);
        }
    }

    valType lookup(const keyType k)
    {
        if(k == key)
        {
            return value;
        }
        else if(k < key)
        {
            return lookupChild(k, left);
        }
        else
        {
            return lookupChild(k, right);
        }
    }

    friend std::ostream &operator<<(std::ostream &out, const BinaryTree &B)
    {
        out << "([" << B.key << "," << B.value << "], " << B.left << ", " << B.right << ")";
        return out;
    }

   friend std::ostream &operator<<(std::ostream &out, const std::unique_ptr<BinaryTree> &Bptr)
   {
       if (Bptr)
       {
           out << (*Bptr);
       }
       else
       {
           out << "(Leaf)";
       }
       return out;
   }

    private:
    /**
     * @brief If child already exists, call recursive insert on child sub-tree. Otherwise create new child node for k-v pair.
     * 
     * @param k key
     * @param v value
     * @param c pointer to child
     */
    void insertChild(const keyType k, const valType v, std::unique_ptr<BinaryTree> &c)
    {
        if(c)
        {
            c->insert(k, v);
        }
        else
        {
            c = std::make_unique<BinaryTree>(k, v);
            c->parent = this;
        }
    }

    /**
     * @brief If child exists, call recursive lookup on child sub-tree. Otherwise throw exception.
     * 
     * @param k key
     * @param c pointer to child node
     * @return valType 
     */
    valType lookupChild(const keyType k, std::unique_ptr<BinaryTree> &c)
    {
        if(c)
        {
            return c->lookup(k);
        }
        else
        {
            throw std::runtime_error("Key value not found.");
        }
    }

    keyType key;
    valType value;

    std::unique_ptr<BinaryTree> left;
    std::unique_ptr<BinaryTree> right;
    BinaryTree *parent;
};