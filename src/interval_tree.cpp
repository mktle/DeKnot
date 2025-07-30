#include "interval_tree.h"

#include <stddef.h>

IntervalNode::IntervalNode(int l, int u){
    this->lower = l;
    this->upper = u;
    this->max = u;
    this->left = NULL;
    this->right = NULL;
}

bool IntervalNode::overlaps(int comp_l, int comp_u){
    return comp_l <= this->upper && comp_u >= this->lower;
}

IntervalNode* IntervalTree::insert(IntervalNode* root, int lower_val, int upper_val){
    if (root == NULL){
        return new IntervalNode(lower_val, upper_val);
    }

    int l = root->lower;

    if (lower_val < l){
        root->left = insert(root->left, lower_val, upper_val);
    }else{
        root->right = insert(root->right, lower_val, upper_val);
    }

    if (root->max < upper_val)
        root->max = upper_val;
    return root;
}

bool IntervalTree::overlap_search(IntervalNode* root, int search_l, int search_u){
    if (root == NULL) return false;
    if (root->overlaps(search_l, search_u)) return true;

    if (root->left != NULL && root->left->max >= search_l){
        return overlap_search(root->left, search_l, search_u);
    }
    return overlap_search(root->right, search_l, search_u);
}

void IntervalTree::destroy_tree(IntervalNode* root){
    if (root == NULL) return;
    destroy_tree(root->left);
    destroy_tree(root->right);
    delete root;
}
