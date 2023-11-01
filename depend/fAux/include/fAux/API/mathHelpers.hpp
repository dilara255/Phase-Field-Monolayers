#pragma once

#include "miscStdHeaders.h"

#define AZ_MAX_UIPOW 1000

namespace AZ{
    template <typename T>
    constexpr T uipow(T num, uint32_t pow)
    {
        return (pow > AZ_MAX_UIPOW ? 0 :
            (pow == 0 ? 1 : num * ipow(num, pow-1));
    }
}
