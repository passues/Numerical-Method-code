import pandas as pd

def quicksort(nums, start, end):
    if start >= end:
        return

    else:
        pivot_index = (start + end)//2
        pivot_value = nums[pivot_index]
        nums[end], nums[pivot_index] = nums[pivot_index], nums[end]
        p = start
        for index in range(start, end):
            if nums[index] < pivot_value:
                nums[index], nums[p] = nums[p], nums[index]
                p += 1
        nums[p], nums[end] = nums[end], nums[p]

        quicksort(nums, start, p -1 ) 
        quicksort(nums, p+1, end)


if __name__=='__main__':
    print("hello world")
    nums = [12,3,21,321,3,12,3,21,4,43,43,5,43,5,345,43,53]
    quicksort(nums, 0, len(nums)-1)
    print(nums)
