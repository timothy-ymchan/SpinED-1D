def digits_str_to_num(digits,base):
    num = 0
    for d in digits:
        intd = int(d)
        assert intd < base, "Contains a digit that is larger than base"
        num = num*base + intd
    return num

def to_base(num, base,ndigit=-1):
    digits = num_to_digits(num,base)
    digits_str = ''.join(str(x) for x in digits[::-1])
    if ndigit > len(digits):
        digits_str = '0'*(ndigit-len(digits)) + digits_str
    return digits_str

def num_to_digits(num,base):
    if num == 0:
        return [0]
    digits = []
    while num:
        digits.append(int(num % base))
        num //= base
    return digits

def circ_shift(num,shift,base,ndit):
    maxval = base**ndit
    assert num < maxval, f"Overflow error: Number does not fit into a {ndit} digit number of base {base}"
    shift = shift % ndit
    head = num // base**(ndit-shift)
    tail = num * base**shift 
    return (head + tail) % maxval

def change_digit(num,d,pos,base,ndigit):
    pos = pos % ndigit
    digit_pos = base**pos
    assert d < base, f"Overflowe Error: Replace digit {d} is larger than base {base}"
    d0 = (num // digit_pos) % base 
    return num + (d-d0)*digit_pos


# def change_digits(original_val,replace_val,di,df,base,ndigit):
#     assert df > di, "The final location must be larger than the initial location"
#     maxval = base**(df-di)
#     assert maxval > replace_val, f"Overflow Error: replace_val={replace_val} exceeds base**(df-di)={maxval} with base={base}, df={df}, di={di}"
#     shift_val = circ_shift(original_val,-di,base,ndigit)
#     original_digits = shift_val % maxval
#     new_val += replace_val - original_digits
#     return circ_shift(new_val,di,base,ndigit)


if __name__ == "__main__":
    num = 11
    print('to_base: ')
    print(to_base(num,base=3,ndigit=5))

    print('circ_shift test:')
    for i in range(10):
        print(to_base(circ_shift(num,i,base=3,ndit=6),base=3,ndigit=6))

    print("digit_str_to_num test")
    digit_str = "2222"
    num = digits_str_to_num(digit_str,3)
    print(f"The number for {digit_str} in base 3 is ",num)
    print(f"The base 3 value of {num} is ",to_base(num,3))

    print('change_digit test:')
    digit_str = "1020122011"
    print(digit_str)
    num = digits_str_to_num(digit_str,3)
    num = change_digit(num,d=2,pos=0,base=3,ndigit=len(digit_str))
    print(to_base(num,3))
    num = change_digit(num,d=0,pos=-3,base=3,ndigit=len(digit_str))
    print(to_base(num,3))