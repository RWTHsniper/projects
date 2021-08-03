#include <iostream>
#include <vector>

// g++ vect.cpp -std=c++11 -o v

class Cat{
public:
    explicit Cat (int age): mAge{age}{}
    void speak() const{
        std::cout << "meow~" << mAge << std::endl;
    }
private:
    int mAge;
};

class Cat2{
public:
    explicit Cat2 (std::string name, int age): mName{std::move(name)},mAge{age}{}
    void speak() const{
        std::cout << "meow~" << mName << " " << mAge << std::endl;
    }
private:
    std::string mName;
    int mAge;
};

class Cat4{
public:
    explicit Cat4 (std::string name): mName{std::move(name)}{
        std::cout << mName << "Cat constructor" << std::endl;
    }
    ~Cat4() noexcept // no exception
    {
        std::cout << mName << "~Cat()" << std::endl;
    }
    Cat4(const Cat4& other): mName(other.mName)
    {
        // compiler tries to copy instead of move by default to avoid exceptions
        std::cout << mName << "copy constructor" << std::endl;
    }
    Cat4(Cat4&& other) noexcept: mName{std::move(other.mName)}
    {
        std::cout << mName << "move Constructor" << std::endl;
    }
private:
    std::string mName;
};

int main()
{
    // vector: container for dynamic array
    // when nums disappears in stack, automatically heap mem is freed
    // std::vector<int> nums(5); // heap memory
    std::vector<int> nums{1,2,3,4};

    nums.emplace_back(5);
    std::cout << nums.size() << std::endl;
    nums.pop_back(); // reduce size
    std::cout << nums.size() << std::endl;

    for(std::size_t idx=0; idx<nums.size(); idx++){
        std::cout << nums[idx] << std::endl;
    }
    for(auto itr=nums.begin(); itr != nums.end(); itr++){
        std::cout << (*itr) << std::endl;
    }

    // ranged for
    for(const int & num: nums){
        std::cout << num << std::endl;
    }

    std::vector<Cat> cats;
    cats.emplace_back(Cat(1));
    cats.emplace_back(Cat(1));
    cats.emplace_back(Cat(1));
    cats.emplace_back(Cat(1));
    cats.emplace_back(Cat(1));

    for(const auto & cat : cats){
        cat.speak();
    }


    std::vector<int> nums2(10000,1);

    nums.emplace_back(2);
    nums.pop_back();

    // vector
    // random access O(1)
    // insertion or removal of elements at the end O(1)
    std::vector<int> nums3(10000,1);
    nums3.emplace_back(2);
    nums3.pop_back();
    // linear in the distance to the end of the vector O(n)
    nums3.emplace(nums3.begin(),3); // to put 3 at the beginning, move every element
    nums3.erase(nums3.begin());

    std::vector<Cat2> cats2;
    cats2.emplace_back(Cat2{"cat0",0}); // unnecessary move operations
    cats2.emplace_back(Cat2{"cat1",1});
    cats2.emplace_back("kitty",2); // just pass arguments for the constructor to avoid declaring temporary variable. Directly declare an object in the vector.
    // Cat& cat = cats2.emplace_back("kitty2",2); // c++17

    for (const auto& cat : cats2){
        cat.speak();
    }

    /* lecture 4
    reserve: reserve capacity to avoid having additional cost to expand a vector
    emplace_back costs move and copy
    */
    std::vector<Cat4> cats4;
    // cats4.reserve(2); // The best way
    cats4.emplace_back("Kitty");
    cats4.emplace_back("nabi"); // Kitty is moved because nabi is inserted

    /* lecture 5
    vector for loop
    How to check performance?
    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "idx Loop:" << diff.count() << "s\n"
    */

    std::cout << "how many elements?" << std::endl;
    std::size_t n;
    std::cin >> n;

    std::vector<int> numsA(n,1);
    std::vector<int> numsB(n,1);
    std::vector<int> numsC(n,1);

    //index based
    for(std::size_t idx=0; idx<numsA.size(); idx++)
    {
        numsA[idx] *= 2;
    }

    //iterator based
    for(auto itr = numsB.begin(); itr != numsB.end(); itr++)
    {
        (*itr) *= 2;
    }

    //range based loop
    for(auto & num : numsC)
    {
        num *= 2;
    }

    // The case where only index based should be used
    // If emplace_back is called, the entire vector moves for extension.
    // If the size of a vector is changing, index-based should be used
    std::vector<int> nums4{0,1,0,1};

    //index based
    for(std::size_t idx=0; idx<nums4.size(); idx++)
    {
        if (nums4[idx]==0)
        {
            nums4.emplace_back(2);
        }
    }

    //iterator based
    for(auto itr = nums4.begin(); itr != nums4.end(); itr++)
    {
        if( (*itr) == 0)
        {
            nums4.emplace_back(2);
        }
    }

    //range based loop
    for(auto & num : nums4)
    {
        if (num == 0)
        {
            nums4.emplace_back(2);
        }
    }

    /*
    Lecture 5
    */


    return 0;

}