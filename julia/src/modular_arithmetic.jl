# Convert a string of digits to a number in a given base
function digits_str_to_num(digits::String, base::Int)
    num::Int = 0
    for ch in digits
        intd = parse(Int, ch)
        @assert intd < base "Contains a digit that is larger than base"
        num = num * base + intd
    end
    return num
end


# Convert a number to a list of digits (LSB first)
function num_to_digits(num::Int, base::Int)
    if num == 0
        return [0]
    end
    digits = Int[]
    while num != 0
        push!(digits, mod(num,base))
        num = div(num, base)
    end
    return digits
end


# Convert a number to a base representation string
function to_base(num::Int, base::Int, ndigit::Int = -1)
    digits = num_to_digits(num, base)
    digits_str = join(reverse(digits))
    if ndigit > length(digits)
        digits_str = repeat("0", ndigit - length(digits)) * digits_str
    end
    return digits_str
end


# Circular digit shift in given base
function circ_shift(num::Int, shift::Int, base::Int, ndit::Int)
    maxval::Int = base^ndit
    @assert num < maxval "Overflow error: Number $num does not fit into a $ndit digit number of base $base (maxval:$maxval)"

    shift = mod(shift,ndit)
    head = div(num, base^(ndit - shift))
    tail = num * base^shift

    return mod((head + tail),maxval)
end


# Replace the digit at index `pos`
function change_digit(num::Int, d::Int, pos::Int, base::Int, ndigit::Int)
    pos = mod(pos,ndigit)
    digit_pos = base^pos

    @assert d < base "Overflow Error: Replace digit $d is larger than base $base"

    d0 = mod((div(num, digit_pos)), base)
    return num + (d - d0) * digit_pos
end

# function main()
#     num = 11
#     println("to_base: ")
#     println(to_base(num, 3, 5))

#     println("\ncirc_shift test:")
#     for i in 0:9
#         println(to_base(circ_shift(num, i, 3, 6), 3, 6))
#     end

#     println("\ndigits_str_to_num test")
#     digit_str = "2222"
#     num = digits_str_to_num(digit_str, 3)
#     println("The number for $digit_str in base 3 is $num")
#     println("The base 3 value of $num is ", to_base(num, 3))

#     println("\nchange_digit test:")
#     digit_str = "1020122011"
#     println(digit_str)
#     num = digits_str_to_num(digit_str, 3)
#     nd = length(digit_str)

#     num = change_digit(num, 2, 0, 3, nd)
#     println(to_base(num, 3))

#     num = change_digit(num, 0, -3, 3, nd)
#     println(to_base(num, 3))
# end

