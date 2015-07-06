#include <vector>
#include <set>

#include <iterator>		// std::back_inserter
#include <algorithm>		// std::set_difference

#include "WireCellUtil/Testing.h"


using namespace WireCell;
using namespace std;
	

int main()
{
    // http://stackoverflow.com/a/10604500
    vector<int> items = {0,1,2,3,4,5,6,7,8,9};
    set<int> die = {2,4,5,42};
    set<int> die2;

    vector<int> result;
    set_difference(items.begin(), items.end(),
		   die.begin(), die.end(),
		   back_inserter(result));

    vector<int> result2;
    set_difference(items.begin(), items.end(),
		   die2.begin(), die2.end(),
		   back_inserter(result2));

    Assert(items.size() == 10);
    Assert(die.size() == 4);
    Assert(result.size() == 7);

    Assert(die2.size() == 0);
    Assert(result2.size() == 10);

}
