#include "./train_kedf.h"

int main()
{
    torch::set_default_dtype(caffe2::TypeMeta::fromScalarType(torch::kDouble));
    auto output = torch::get_default_dtype();
    std::cout << "Default type: " << output << std::endl;

    Train_KEDF train;
    train.input.readInput();
    if (train.input.check_pot)
    {
        train.potTest();
    }
    else
    {
        train.init();
        train.train();
    }
}