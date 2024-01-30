"""
    index_generator()

Create an index generator that yields unique 3D indices in a sorted order based on the Euclidean distance 
from the origin. The generator maintains a queue of indices and dynamically calculates new indices as needed.

# Returns
- `yield_index()`: A function to retrieve the next unique 3D index.
- `reset()`: A function to reset the generator by emptying the index queue.

## Example

```julia
# Create an index generator
gen_yield, gen_reset = index_generator()

# Get the next index
index1 = gen_yield()  # e.g., (0, 0, 1)
index2 = gen_yeild()  # e.g., (0, 1, 0)

# Reset the generator
gen_reset()

# Get another index after reset
index3 = gen_yield()  # e.g., (0, 0, 1)
"""
function index_generator()
    index_queue = Queue{Tuple{Int32, Int32, Int32}}()
    max_index = 0
    function yield_index()
        if length(index_queue) == 0
            max_index += 1
            calculate_indices(max_index)
        end
        dequeue!(index_queue)
    end
    function calculate_indices(max_index)
        indices = vcat(
                collect(
                    Iterators.product(
                        -max_index:max_index, 
                        -max_index:max_index, 
                        -max_index:max_index)
                        )
                    ...)
        
        sorted_indices = sort(indices, by=x->sum(x.^2))[2:end]
        enqueue!.(Ref(index_queue), sorted_indices)
    end
    function reset()
        empty!(index_queue)
    end
    return yield_index, reset
end