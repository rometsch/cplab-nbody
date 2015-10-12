/*
 * nbodysys_test.cpp
 *
 *  Created on: Oct 12, 2015
 *      Author: thomas
 */

#include "nbsys.h"

#include "gtest/gtest.h"


TEST(LoadTest, Correctness) {
	nbsys *s = new nbsys("./data/pla3");
	EXPECT_EQ(3,s->N);
	EXPECT_DOUBLE_EQ(220,s->t_max);
	EXPECT_DOUBLE_EQ(0.1,s->eta);
	EXPECT_DOUBLE_EQ(0.9999,s->m[0]);
	EXPECT_DOUBLE_EQ(0.00001,s->m[1]);
	EXPECT_DOUBLE_EQ(0.00009,s->m[2]);
	EXPECT_DOUBLE_EQ(-0.0001925,s->Rx);
	EXPECT_DOUBLE_EQ(4.9999994000000006e-05,s->Py);
	EXPECT_DOUBLE_EQ(0,s->Pz);
	EXPECT_DOUBLE_EQ(-1-4.9999994000000006e-05,s->vy[1]);
	EXPECT_DOUBLE_EQ(-2.25+0.0001925,s->rx[2]);
}
