#ifndef ATEN_KERNELS_EINSUM_OP_H_
#define ATEN_KERNELS_EINSUM_OP_H_

#include <ATen/core/tensor.h>
#include <ATen/core/tensor_map.h>

namespace container {

struct EinsumOption {
    bool conj_x = false;
    bool conj_y = false;
    float alpha = 1.0;
    float beta  = 0.0;
    Tensor* out = nullptr;

    EinsumOption(bool conj_x_ = false, bool conj_y_ = false, float alpha_ = 1.0, float beta_ = 0.0, Tensor* out_ = nullptr)
        : conj_x(conj_x_), conj_y(conj_y_), alpha(alpha_), beta(beta_), out(out_) {}
};

namespace einsum_utils {
struct BCast;

// Dummy axis label used to denote an ellipsis in an input or output subscript.
constexpr int kEllipsisLabel = -1;

// Each dimension is categorized into exactly one of five types based on
// whether its corresponding label is present in the input and/or the output
// subscripts.
enum EinsumDimensionType {
    // Batch dimensions are those present in two inputs as well as the output.
    // They are part of the batch dimensions during Tensor contraction. Such
    // dimensions may be broadcasting dimensions (those mapping to ellipsis)
    // or explicit batch dimensions corresponding to named axis labels.
    kBroadcasting = 0,
    kBatch = 1,
    // Free dimensions are present in exactly one of the inputs, and also the
    // output. These are non-contracted axes in the Tensor contraction.
    kFree = 2,
    // Contract dimensions are present in two inputs, but not the output. These
    // dimensions are contracted in Tensor contraction.
    kContract = 3,
    // Reduce dimensions are present in exactly one input; and not in the output
    // and are summed over prior to Tensor contraction.
    kReduce = 4,
};

// Parses and validates an einsum equation in explicit form.
bool ValidateEinsumEquation(
    const std::string& equation,
    std::vector<std::string>& input_subscripts,
    std::string& output_subscript);

// Parses and validates the equation and the input shapes. Single character
// labels are integerized, and we populate input and output label subscripts
// and corresponding counts. Also create the mapping from (named) labels to
// their EinsumDimensionType.
bool ParseEinsumEquation(
    const std::string& equation, 
    std::vector<EinsumDimensionType>& label_types,
    std::vector<std::vector<int>>& input_labels,
    std::vector<int>& output_labels,
    std::vector<std::vector<int>>& input_label_counts,
    std::vector<int>& output_label_counts,
    std::vector<bool>& input_has_ellipsis,
    bool& output_has_ellipsis);

bool ProcessDimensions(
    const std::vector<const Tensor*>& inputs,
    std::vector<EinsumDimensionType>& label_types,
    std::vector<std::vector<int>>& input_labels,
    std::vector<int>& output_labels,
    std::vector<std::vector<int>>& input_label_counts,
    std::vector<int>& output_label_counts,
    const std::vector<bool>& input_has_ellipsis,
    const bool output_has_ellipsis,
    std::unordered_map<int, int64_t>& label_to_dim_sizes);

// This function records the mapping of a label to its corresponding dimension for a specific axis in the input tensor.
// It also validates that the label and dimension mapping is consistent with previous recordings, ensuring that the 
// same label is not mapped to different dimensions along different axes.
bool RecordLabelToDimension(
    const int label,
    const int axis,
    const Tensor& input,
    std::unordered_map<int, int64_t>& label_to_dim_sizes);

bool ReduceOperand(
    const Tensor& input,
    const std::vector<EinsumDimensionType>& label_types,
    std::vector<int>& labels,
    const std::vector<int>& label_counts,
    std::vector<int>& free_labels,
    int& swap_free_and_contract,
    Tensor& output);

/**
 * @brief A function to perform contraction operation on multiple Tensors.
 *
 * This functor applies a contraction operation on multiple input Tensors and computes the result.
 * The contraction operation combines the input Tensors based on a specific contraction pattern to produce
 * a single output Tensor. The contraction pattern is defined by the `swap_free_and_contract` vector, which
 * specifies whether each input Tensor should be contracted or simply copied to the output.
 *
 */
bool ContractOperands(
    std::vector<Tensor>& inputs,
    const std::vector<int>& swap_free_and_contract,
    const EinsumOption& option,
    Tensor& output);

void ProcessOutput(
    const Tensor& input,
    const std::vector<einsum_utils::EinsumDimensionType>& label_types,
    const std::vector<std::vector<int>>& free_labels,
    std::unordered_map<int, int64_t>& label_to_dim_sizes,
    const std::vector<int>& output_labels,
    const std::vector<int>& output_label_counts,
    Tensor& output);

} // namespace einsum_utils

namespace op {

// TODO: implement this method this week!
// piapia pat face

/**
 *
 * @brief Computes the Einstein summation convention on the given tensors
 * based on the expression passed as a string.
 *
 * @tparam Tensors Variadic template parameter pack for the input tensors.
 *
 * @param equation The Einstein summation convention expression.
 * @param tensors The input tensors for the summation operation.
 * @return The resulting tensor after applying the Einstein summation convention on the input tensors.
 *
 * @throws std::invalid_argument if the expression or the input tensors are invalid.
 * @throws std::runtime_error if an error occurs while performing the summation operation.
 */
template <typename... Tensors>
typename std::enable_if<std::is_same<
        typename std::common_type<Tensors...>::type, Tensor>::value, Tensor>::type
    einsum_impl(const std::string& equation, const EinsumOption& option, const Tensors&... tensors)
{
    // Check the input dimension
    constexpr int num_inputs = sizeof...(Tensors);
    if (num_inputs > 2) {
        throw std::invalid_argument("Einstein notation only support two or less tensors!");
    }
    const std::vector<const Tensor*> inputs{reinterpret_cast<const Tensor*>(&tensors)...};
    // Init the input and output labels
    std::vector<std::vector<int>> input_labels = {};
    std::vector<int> output_labels = {};
    std::vector<einsum_utils::EinsumDimensionType> label_types = {};
    std::vector<std::vector<int>> input_label_counts = {};
    std::vector<int> output_label_counts = {};
    std::vector<bool> input_has_ellipsis = {};
    bool output_has_ellipsis = {};

    einsum_utils::ParseEinsumEquation(
        equation, label_types, 
        input_labels, output_labels, 
        input_label_counts, output_label_counts, 
        input_has_ellipsis, output_has_ellipsis);
    
    if (input_labels.size() != num_inputs) {
        throw std::runtime_error("The number of input tensors does not match the number of input labels!");
    }
    
    std::unordered_map<int, int64_t> label_to_dim_sizes = {};

    einsum_utils::ProcessDimensions(
        inputs, label_types,
        input_labels, output_labels, 
        input_label_counts, output_label_counts,
        input_has_ellipsis, output_has_ellipsis,
        label_to_dim_sizes);

    std::vector<std::vector<int>> free_labels(num_inputs);
    std::vector<int> swap_free_and_contract(num_inputs);
    std::vector<Tensor> inputs_reduced(num_inputs, Tensor(DataType::DT_FLOAT, {}));

    for (int ii = 0; ii < num_inputs; ++ii) {
        einsum_utils::ReduceOperand(
            *inputs[ii], label_types,
            input_labels[ii], input_label_counts[ii],
            free_labels[ii], swap_free_and_contract[ii], 
            inputs_reduced[ii]);
    }

    // After reduction, the inputs should be reshaped to Tensors suitable for
    // contraction. If num_inputs is 1, the reduced input is simply forwarded to
    // the output.
    Tensor contraction_output_reshaped;
    einsum_utils::ContractOperands(
        inputs_reduced, swap_free_and_contract, 
        option, contraction_output_reshaped);
    
    Tensor output;
    // Copy the batch labels from the contraction output. Recover the batch
    // shape, which may have been broadcasted.
    einsum_utils::ProcessOutput(
        contraction_output_reshaped, label_types,
        free_labels, label_to_dim_sizes, 
        output_labels, output_label_counts,
        output);

    return std::move(output);
}

// Make the conj params only works for the matmul equations.
inline Tensor einsum(const std::string& equation, const Tensor& A) {
    const EinsumOption& option = {};
    return std::move(op::einsum_impl(equation, option, A));
}

inline Tensor einsum(const std::string& equation, const Tensor& A, const Tensor& B, const EinsumOption& option = {}) {
    return std::move(op::einsum_impl(equation, option, A, B));
}

} // namespace op
} // namespace container

#endif // ATEN_KERNELS_EINSUM_OP_H_