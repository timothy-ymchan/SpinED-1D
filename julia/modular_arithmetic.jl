
function circ_shift(num::Int,shift::Int;base::Int,ndit::Int)
    maxval = base^ndit 
    num < maxval || throw(DomainError(num,"Number does not fit into a $ndit digit number of base $base"))
    shift = mod(shift,ndit)
    cut = base^(ndit-shift)
    head, tail = divrem(num,cut)
    return mod(head + tail*base^shift,maxval)
end

function change_digit(num::Int,digit::Int,pos::Int;base::Int,ndit::Int)
    pos = base^mod(pos,ndit)
    BoundsError
    (0<= digit && digit< base) || throw(DomainError(digit,"Replace digit $digit is larger than base $base"))
    d0 = mod(div(num,pos),base)
    return num + (digit-d0)*pos
end

function digits_to_num(digits::String,base::Int)
    return parse(Int,digits;base=base)
end

function num_to_digits(num::Int,base::Int;ndigit::Int=-1)
    return lpad(string(num,base=base),ndigit,'0') # Pad with zeros if necessary
end


# let 
#     base=3
#     num = 11 
#     @show num_to_digits(num,3;ndigit=5)
#     println("circ_shift test")
#     for i in 0:10
#         println(num_to_digits(circ_shift(num,i;base=base,ndit=6),base,ndigit=6))
#     end

#     println("digits_to_num test")
#     digit = "2222"
#     num = digits_to_num(digit,base)
#     println("The number for $digit in base 3 is $num")
#     println("The base 3 digits for $num is $(num_to_digits(num,base))")

#     println("Change digit test")
#     digit = "1020122011"
#     println(digit)
#     num = digits_to_num(digit,base)
#     num = change_digit(num,2,0;base=3,ndit=length(digit))
#     println(num_to_digits(num,base))
#     num = change_digit(num,0,-3;base=3,ndit=length(digit))
#     println(num_to_digits(num,base))
# end