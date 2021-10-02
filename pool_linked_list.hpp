#ifndef ROLLER_POOL_LINKED_LIST_HPP_
#define ROLLER_POOL_LINKED_LIST_HPP_

#include <memory>

namespace roller {

/**
  Please don't insert more than `maxSize` elements. It might segfault.
  Also don't remove the head. It probably won't segfault, but it will misbehave.
*/
template <class T> class PoolLinkedList {

public:
  struct Node {

    T obj;
    Node *prev;
    Node *next;
  };

private:
  std::unique_ptr<Node[]> pool;
  std::unique_ptr<std::size_t[]> stack;
  std::size_t ind, maxSize;
  Node *head;

  Node *fetchFromPool() {
    std::size_t i = stack[--ind];
    return pool.get() + i;
  }

  void returnToPool(Node *node) {
    std::size_t i = node - pool.get();
    stack[ind++] = i;
  }

public:
  explicit PoolLinkedList(std::size_t maxSize)
      : pool(new Node[maxSize + 1]), stack(new std::size_t[maxSize]),
        ind(maxSize), maxSize(maxSize), head(pool.get()) {
    for (std::size_t i = maxSize, j = 0; i--; j++) {
      stack[j] = i + 1;
    }
    head->prev = head;
    head->next = head;
  }

  PoolLinkedList()
      : pool(new Node[9]), stack(new std::size_t[8]), ind(8), head(pool.get()) {
    for (std::size_t i = 8, j = 0; i--; j++) {
      stack[j] = i + 1;
    }
    head->prev = head;
    head->next = head;
  }

  void clear() {
    ind = maxSize;
    for (std::size_t i = maxSize, j = 0; i--; j++) {
      stack[j] = i + 1;
    }
    head->prev = head;
    head->next = head;
  }

  struct iterator {

    Node *n;

    typedef std::ptrdiff_t difference_type;
    typedef T value_type;
    typedef T *pointer;
    typedef T &reference;
    typedef std::bidirectional_iterator_tag iterator_category;

    reference operator*() const noexcept { return n->obj; }
    pointer operator->() const noexcept { return &n->obj; }
    iterator &operator++() noexcept {
      n = n->next;
      return *this;
    }
    iterator operator++(int) noexcept {
      iterator toreturn = *this;
      operator++();
      return toreturn;
    }
    iterator &operator--() noexcept {
      n = n->prev;
      return *this;
    }
    iterator operator--(int) noexcept {
      iterator toreturn = *this;
      operator--();
      return toreturn;
    }

    bool operator==(const iterator &other) const noexcept {
      return n == other.n;
    }
    bool operator!=(const iterator &other) const noexcept {
      return n != other.n;
    }
  };

  typedef const iterator const_iterator;

  iterator begin() { return {head->next}; }
  const_iterator begin() const { return {head->next}; }
  const_iterator cbegin() const { return {head->next}; }
  iterator end() { return {head}; }
  const_iterator end() const { return {head}; }
  const_iterator cend() const { return {head}; }
  iterator rbegin() { return {head->prev}; }
  const_iterator rbegin() const { return {head->prev}; }
  const_iterator crbegin() const { return {head->prev}; }
  iterator rend() { return {head}; }
  const_iterator rend() const { return {head}; }
  const_iterator crend() const { return {head}; }

  typedef T value_type;
  typedef T &reference;
  typedef const T &const_reference;
  typedef typename iterator::difference_type difference_type;
  typedef std::size_t size_type;

  Node *getHead() { return head; }
  const Node *getHead() const { return head; }
  Node *front() { return head->next; }
  const Node *front() const { return head->next; }
  Node *back() { return head->prev; }
  const Node *back() const { return head->prev; }

  std::size_t size() const { return maxSize - ind; }

  std::size_t invSize() const { return ind; }

  Node *insertBefore(Node *node, const T &obj) {
    Node *nn = fetchFromPool();
    Node *p = node->prev;
    nn->obj = obj;
    nn->prev = p;
    nn->next = node;
    p->next = nn;
    node->prev = nn;
    return nn;
  }
  Node *insertBefore(Node *node, T &&obj) {
    Node *nn = fetchFromPool();
    Node *p = node->prev;
    nn->obj = obj;
    nn->prev = p;
    nn->next = node;
    p->next = nn;
    node->prev = nn;
    return nn;
  }
  Node *insertAfter(Node *node, const T &obj) {
    Node *nn = fetchFromPool();
    Node *n = node->next;
    nn->obj = obj;
    nn->prev = node;
    nn->next = n;
    n->prev = nn;
    node->next = nn;
    return nn;
  }
  Node *insertAfter(Node *node, T &&obj) {
    Node *nn = fetchFromPool();
    Node *n = node->next;
    nn->obj = obj;
    nn->prev = node;
    nn->next = n;
    n->prev = nn;
    node->next = nn;
    return nn;
  }

  Node *addToBeg(const T &obj) { return insertAfter(head, obj); }
  Node *addToBeg(T &&obj) { return insertAfter(head, obj); }
  Node *addToEnd(const T &obj) { return insertBefore(head, obj); }
  Node *addToEnd(T &&obj) { return insertBefore(head, obj); }

  void remove(Node *node) {
    Node *p = node->prev;
    Node *n = node->next;
    p->next = n;
    n->prev = p;
    returnToPool(node);
  }

  void removeFromBeg() { remove(head->next); }
  void removeFromEnd() { remove(head->prev); }
};

} // namespace roller

#endif // ROLLER_POOL_LINKED_LIST_HPP_

