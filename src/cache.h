#pragma once

#include <list>
#include <unordered_map>

template<class Key, class Value>
class Cache
{
public:
	typedef Key key_type;
	typedef Value value_type;

	explicit Cache(size_t max_size):maxSize(max_size)
	{
		data.reserve(max_size);
	}
	
	const value_type* find(key_type k)
	{
		auto it = data.find(k);
		if (it == data.end()) return nullptr;
		
		lru.splice(lru.begin(), lru, it->second);
		return &it->second->second;
	}
	
	void insert(key_type k, const value_type& value)
	{
		lru.emplace_front(k, value);
		auto res = data.emplace(k, lru.begin());
		
		if (!res.second)
		{
			lru.erase(res.first->second);
			res.first->second = lru.begin();
		}
		
		while (data.size() > maxSize)
		{
			data.erase(lru.back().first);
			lru.pop_back();
		}
	}
	
private:
	typedef std::list< std::pair< key_type, value_type> > LRU;
	
	size_t maxSize;
	LRU lru;
	std::unordered_map< key_type, typename LRU::iterator > data;
};
