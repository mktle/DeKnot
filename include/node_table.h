#ifndef NODE_TABLE_H
#define NODE_TABLE_H

#include <cstddef>
#include <cstdint>
#include <vector>

class Node;

class NodeTable{
    public:
        NodeTable() : mask(0), count(0), zero_value(NULL), has_zero(false) {}

        size_t size() const { return this->count + (this->has_zero ? 1 : 0); }

        void reserve(size_t n){
            size_t need = 16;
            while (need < 2 * n + 1) need <<= 1;
            if (need > this->slots.size()) this->rehash(need);
        }

        Node* find(uint64_t key) const {
            if (key == 0) return this->has_zero ? this->zero_value : NULL;
            if (this->slots.empty()) return NULL;
            size_t i = this->index_of(key);
            while (true){
                const Slot& s = this->slots[i];
                if (s.key == key) return s.value;
                if (s.key == 0) return NULL;
                i = (i + 1) & this->mask;
            }
        }

        void insert_or_assign(uint64_t key, Node* value){
            if (key == 0){
                this->zero_value = value;
                this->has_zero = true;
                return;
            }
            if ((this->count + 1) * 2 > this->slots.size()) this->rehash(this->slots.empty() ? 16 : this->slots.size() * 2);
            size_t i = this->index_of(key);
            while (true){
                Slot& s = this->slots[i];
                if (s.key == key){
                    s.value = value;
                    return;
                }
                if (s.key == 0){
                    s.key = key;
                    s.value = value;
                    this->count++;
                    return;
                }
                i = (i + 1) & this->mask;
            }
        }

        void erase(uint64_t key){
            if (key == 0){
                this->has_zero = false;
                this->zero_value = NULL;
                return;
            }
            if (this->slots.empty()) return;
            size_t i = this->index_of(key);
            while (true){
                Slot& s = this->slots[i];
                if (s.key == 0) return;
                if (s.key == key) break;
                i = (i + 1) & this->mask;
            }
            // backward shift
            size_t hole = i;
            size_t j = (i + 1) & this->mask;
            while (this->slots[j].key != 0){
                size_t home = this->index_of(this->slots[j].key);
                // entry j moves to the hole only if the hole lies between its home and j
                bool movable = (hole <= j) ? (home <= hole || home > j) : (home <= hole && home > j);
                if (movable){
                    this->slots[hole] = this->slots[j];
                    hole = j;
                }
                j = (j + 1) & this->mask;
            }
            this->slots[hole].key = 0;
            this->slots[hole].value = NULL;
            this->count--;
        }

    private:
        struct Slot {
            uint64_t key;
            Node* value;
        };
        std::vector<Slot> slots;
        size_t mask;
        size_t count;
        Node* zero_value;
        bool has_zero;

        size_t index_of(uint64_t key) const {
            uint64_t x = key;
            x ^= x >> 30;
            x *= 0xbf58476d1ce4e5b9ULL;
            x ^= x >> 27;
            x *= 0x94d049bb133111ebULL;
            x ^= x >> 31;
            return (size_t)x & this->mask;
        }

        void rehash(size_t new_size){
            std::vector<Slot> old;
            old.swap(this->slots);
            this->slots.assign(new_size, Slot{0, NULL});
            this->mask = new_size - 1;
            this->count = 0;
            for (size_t i{0}; i < old.size(); i++){
                if (old[i].key != 0) this->insert_or_assign(old[i].key, old[i].value);
            }
        }
};

#endif
