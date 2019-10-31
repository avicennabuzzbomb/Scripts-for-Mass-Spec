### playing around with decimal rounding!
if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it
    import numpy
    import decimal

    val1 = decimal.Decimal(632.4556001) / decimal.Decimal(1)   # when typecasting to Decimal, need to create a 'Decimal' instance of Class decimal.
    val2 = decimal.Decimal(632.456) / decimal.Decimal(1)       # when typecasting to Decimal, need to create a 'Decimal' instance of Class decimal.
    val3 = 2.6655
    val4 = 2.6755

    """
    -Example usage-
    a = decimal.Decimal('5.1')
    b = decimal.Decimal('3.14')
    c = 4
    d = 3.14
    """

    #before rounding
    print("Before round: ", val1, " ", val2, " ", val3, " ", val4)

    #after rounding
    new_val1 = round(val1, 3)
    new_val2 = round(val2, 3)
    new_val3 = round(val3, 3)
    new_val4 = round(val4, 3)
    print("After round: ", new_val1, " ", new_val2, " ", new_val3, " ", new_val4)
    
    #after rounding all as Decimal type
    val3 = decimal.Decimal(val3)
    val4 = decimal.Decimal(val4)

    new_val1 = round(val1, 3)
    new_val2 = round(val2, 3)
    new_val3 = round(val3, 3)
    new_val4 = round(val4, 3)
    print("After round as Decimal: ", new_val1, " ", new_val2, " ", new_val3, " ", new_val4)

    new = 7.325
    new1 = Decimal('7.325').quantize(Decimal('.01'), rounding=ROUND_DOWN)
    print(new,"rounds to",new1)