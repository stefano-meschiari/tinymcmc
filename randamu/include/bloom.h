/*! \file
 *
 * \brief A probabilistic set data structure.
 */

#ifndef RD_BLOOM_H_
#define RD_BLOOM_H_

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Probabilistic set data structure.
 */
typedef struct {
} rd_bloom;

#ifdef __cplusplus
}
#endif

// It has a low probability of returning false positives (it says an element is in while it's not).

#endif
