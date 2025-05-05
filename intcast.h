#pragma once /*
inline int float2int(double d)
{
    union Cast
    {
        double d;
        long l;
    };
    volatile Cast c;
    c.d = d + 6755399441055744.0;
    return c.l;
} */

// this is the same thing but it's
// not always optimizer safe
inline int float2int(double d)
{
    d += 6755399441055744.0;
    return reinterpret_cast<int&>(d);
}