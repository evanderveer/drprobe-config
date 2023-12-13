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
        
        sorted_indices = sort(indices, by=x->sum(x.^2))
        enqueue!.(Ref(index_queue), sorted_indices)
    end
    function reset()
        empty!(index_queue)
        max_index = 0;
    end
    return yield_index, reset
end