class IntervalNode{
    public:
        IntervalNode(int l, int u);
        bool overlaps(int comp_l, int comp_u);

        int lower;
        int upper;
        int max;
        IntervalNode *left, *right;
};

namespace IntervalTree{
    IntervalNode* insert(IntervalNode *root, int lower, int upper);
    bool overlap_search(IntervalNode* root, int search_l, int search_u);
}
