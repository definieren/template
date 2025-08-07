可删堆。

`Heap<T, cmp>`：元素类型为 `T`，比较方式为 `cmp` 的可删堆。

- `reduce()`：将堆中懒惰删除的元素清空。
- `push(x)`：向堆中加入元素 $x$。
- `pop(x)`：弹出元素（若无传参则为弹出堆顶）。
- `top()`：返回堆顶元素。
- `size()`：返回堆中元素个数。
- `empty()`：返回堆是否为空。

```cpp
template<class T, class cmp = std::less<T>>
struct Heap {
private:
	std::priority_queue<T, std::vector<T>, cmp> addQ, delQ;
public:
	Heap(): addQ{}, delQ{} {}
	void reduce() {
		assert(addQ.size() >= delQ.size());
		while (delQ.size() && addQ.top() == delQ.top()) {
			addQ.pop(), delQ.pop();
		}
	}
	void push(T x) {
		addQ.push(x);
	}
	void pop(T x) {
		delQ.push(x);
	}
	void pop() {
		reduce();
		assert(addQ.size());
		addQ.pop();
	}
	T top() {
		reduce();
		assert(addQ.size());
		return addQ.top();
	}
	size_t size() {
		assert(addQ.size() >= delQ.size());
		return addQ.size() - delQ.size();
	}
	bool empty() {
		return !size();
	}
};
```