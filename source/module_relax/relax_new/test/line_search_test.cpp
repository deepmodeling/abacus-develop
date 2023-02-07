#include "gtest/gtest.h"

/************************************************
 *  unit test of class Line_Search
 ***********************************************/

/**
 * - Tested Functions:
 *   - Line_Search::line_search()
 *   - Line_Search::first_order()
 *   - Line_Search::third_order()
 *   - Line_Search::init_brent()
 *   - Line_Search::update_brent()
 *   - Line_Search::brent()
 */


#define private public
#include "module_relax/relax_new/line_search.h"

class LineSearchPrepare
{
public:
	LineSearchPrepare(
			bool restart_in,
			double x_in,
			double y_in,
			double f_in,
			double xnew_in,
			double conv_thr_in):
		restart(restart_in),
		x(x_in),
		y(y_in),
		f(f_in),
		xnew(xnew_in),
		conv_thr(conv_thr_in)
	{}
	bool restart;
	double x;
	double y;
	double f;
	double xnew;
	double conv_thr;
};


class LineSearchTest : public ::testing::TestWithParam<LineSearchPrepare>{};

TEST_P(LineSearchTest,LineSearch)
{
	LineSearchPrepare lsp = GetParam();
	Line_Search ls;
	EXPECT_EQ(ls.ls_step,0);
	bool test_conv = ls.line_search(lsp.restart,lsp.x,lsp.y,lsp.f,lsp.xnew,lsp.conv_thr);
	EXPECT_EQ(ls.ls_step,1);
	while (!test_conv)
	{
		test_conv = ls.line_search(false,lsp.x,lsp.y,lsp.f,lsp.xnew,lsp.conv_thr);
	}
	EXPECT_GE(ls.ls_step,3);
	EXPECT_TRUE(test_conv);
}

INSTANTIATE_TEST_SUITE_P(LineSearch,LineSearchTest,::testing::Values(
			// restart, x, y, f, xnew, conv_thr
			LineSearchPrepare(true,5.0,100.1,-10.0,0.0,1e-6)
			));

#undef private
