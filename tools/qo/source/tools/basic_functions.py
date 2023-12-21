import numpy as np
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
        else:
            if read_real:
                real += letter
            elif read_imag:
                imag += letter
    try:
        real = float(real)
        imag = float(imag)
        result = real+imag*1j
        return result
    except ValueError:
        print(f"Error: {number}")
        raise ValueError

if __name__ == "__main__":

    print(make_complex("(4.32385486762811e-02,-1.59403849290982e-02)"))