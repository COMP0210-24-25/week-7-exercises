#include <functional>
#include <ostream>
#include <stdexcept>
#include <vector>

using std::vector;

/**
 * @brief Table used for mapping key value pairs. 
 * Each row is indexed by a hash function on the key. Key-value pairs which map to the same row are stored as an unsorted list.
 * 
 * @tparam keyType 
 * @tparam valType 
 */
template<typename keyType, typename valType>
class Table
{
    struct kv_pair
    {
        keyType key;
        valType val;
    };

    public:
    Table(const int N) : size(N)
    {
        table = vector< vector<kv_pair> >(size);
    }

    void insert(const keyType key, const valType val)
    {
        kv_pair kv{key, val};
        auto &row = table[index(key)];
        if (row.empty())
        {
            row.push_back(kv);
        }
        else
        {
            for(auto &row_entry : row)
            {
                if(row_entry.key == key)
                {
                    row_entry.val = val;
                    return;
                }
            }
            row.push_back(kv);
        }
    }

    valType lookup(const keyType key)
    {
        const auto &row = table[index(key)];

        for (auto &row_entry : row)
        {
            if (row_entry.key == key)
            {
                return row_entry.val;
            }
        }
        throw std::runtime_error("Key not found.");
    }

    friend std::ostream &operator<<(std::ostream &out, const Table &h)
    {
        for(int i = 0; i < h.size; i++)
        {
            const auto row = h.table[i];
            if(!row.empty())
            {
                out << i << ": ";
                for(const auto entry : row)
                {
                    out << "(" << entry.key << "," << entry.val << ") ";
                }
                out << std::endl;
            }
        }
        return out;
    }

    private:

    /**
     * @brief Constant time hash function modulo size of table
     * 
     * @param key 
     * @return int 
     */
    int index(keyType key)
    {
        return hash(key) % size;
    }

    std::hash<keyType> hash;
    const int size;
    vector< vector<kv_pair>> table;
};