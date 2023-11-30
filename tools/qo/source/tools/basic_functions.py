def make_complex(number: str) -> complex:
    """make complex number from string

    Args:
        number (str): string of number

    Returns:
        complex: complex number
    """

    read_real = False
    read_imag = False
    real = ""
    imag = ""
    for letter in number:
        if letter == "(":
            read_real = True
            continue
        elif letter == ",":
            read_real = False
            read_imag = True
            continue
        elif letter == ")":
            read_imag = False
            continue
        if read_real:
            real += letter
        elif read_imag:
            imag += letter
    try:
        result = complex(float(real), float(imag))
        return result
    except ValueError:
        print(f"Error: {number}")
        raise ValueError
