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

if __name__ == "__main__":
    num = 11
    print('to_base: ')
    print(to_base(num,base=3,ndigit=5))

    print('circ_shift test:')
    for i in range(10):
        print(to_base(circ_shift(num,i,base=3,ndit=6),base=3,ndigit=6))
