#pragma once

#include <afxwin.h>

#include <vector>

#include <util/common/math/complex.h>
#include <util/common/plot/plot.h>
#include <util/common/math/fft.h>
#include <util/common/math/vec.h>
#include <util/common/math/raster.h>
#include <util/common/geom/geom.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // !M_PI


namespace model
{

    /*****************************************************/
    /*                     params                        */
    /*****************************************************/

    struct parameters
    {
        // system params
        double snr;
        size_t n;
    };

    inline parameters make_default_parameters()
    {
        parameters p =
        {
            // system params
            1.0,
            100,
        };
        return p;
    }

    inline plot::drawable::ptr_t make_root_drawable
    (
        std::vector < plot::drawable::ptr_t > layers
    )
    {
        using namespace plot;

        return viewporter::create(
            layer_drawable::create(layers),
            make_viewport_mapper(world_t{-1, 1, -1, 1})
        );
    }

    inline plot::drawable::ptr_t make_bmp_plot(CBitmap & b)
    {
        return plot::custom_drawable::create([&b] (CDC & dc, const plot::viewport & vp)
        {
            if (!b.m_hObject) return;
            CDC memDC; memDC.CreateCompatibleDC(&dc);
            memDC.SelectObject(&b);
            dc.SetStretchBltMode(HALFTONE);
            auto wh = b.GetBitmapDimension();
            dc.StretchBlt(vp.screen.xmin, vp.screen.ymin,
                          vp.screen.width(), vp.screen.height(),
                          &memDC, 0, 0, wh.cx, wh.cy, SRCCOPY);
        });
    }

    /*****************************************************/
    /*                     data                          */
    /*****************************************************/

    struct bitmap
    {
        size_t h, w;
        std::vector < std::vector < double > > data;
        
        bitmap & reshape(size_t n, size_t m)
        {
            h = n; w = m;
            data.clear();
            data.resize(n, std::vector < double > (m, 0));
            return *this;
        }
        bitmap & to_cbitmap(CBitmap & bmp, double d, bool abs)
        {
            std::vector < COLORREF > buf(h * w);
            for (size_t i = 0; i < h; ++i)
            for (size_t j = 0; j < w; ++j)
            {
                double v = abs ? std::abs(data[i][j]) : std::fmax(0, data[i][j]);
                BYTE c = (BYTE) std::fmin(v / d * 255, 255);
                buf[w * i + j] = RGB(c, c, c);
            }
            bmp.DeleteObject();
            bmp.CreateBitmap(w, h, 1, sizeof(COLORREF) * 8, (LPCVOID) buf.data());
            bmp.SetBitmapDimension(w, h);
            return *this;
        }
        bitmap & from_cbitmap(CBitmap & bmp)
        {
            BITMAP bmpd;
            bmp.GetBitmap(&bmpd);
            w = bmpd.bmWidth; h = bmpd.bmHeight;
            reshape(h, w);
            CDC dc; dc.CreateCompatibleDC(nullptr);
            dc.SelectObject(&bmp);
            for (size_t i = 0; i < h; ++i)
            for (size_t j = 0; j < w; ++j)
            {
                auto p = dc.GetPixel(j, i);
                data[i][j] = (0.299 * (p & 0xff) +
                              0.587 * ((p >> 8) & 0xff) +
                              0.114 * ((p >> 16) & 0xff)) / 255.;
            }
            return *this;
        }
    };

    /*****************************************************/
    /*                     algo                          */
    /*****************************************************/

    inline void make_line(
        size_t w, size_t h, int d, double theta,
        geom::point < int > & p1,
        geom::point < int > & p2)
    {
        using arr4_t = std::array < geom::point2d_t, 4 > ;
        const arr4_t points =
        {
            { { 0, 0 }, { w - 1, 0 }, { w - 1, h - 1 }, { 0, h - 1 } }
        };
        auto rect = geom::make_convex_polygon(points);

        math::v3 < double > c(w / 2.0, h / 2.0);
        math::v3 < double > v1(-1.0 * w, h / 2.0 + d);
        math::v3 < double > v2(2.0 * w, h / 2.0 + d);

        math::m3 < double > rm = math::rotate_z(theta);

        v1 = (rm * (v1 - c) + c);
        v2 = (rm * (v2 - c) + c);

        geom::point2d_t w1 = { v1.x, v1.y };
        geom::point2d_t w2 = { v2.x, v2.y };

        auto intersections = rect.intersections(geom::make_line(w1, w2));

        if (intersections.empty()) p1 = p2 = { 0, 0 };
        else
        {
            p1 = { std::round(intersections[0].x), std::round(intersections[0].y) };
            if (intersections.size() == 1) p2 = p1;
            else p2 = { std::round(intersections[1].x), std::round(intersections[1].y) };
        }
    }

    class backprojection
    {
    private:
        const parameters & p;
        std::vector < math::complex < > > Q, P;
        size_t w, h, S, C;
    public:
        backprojection(const parameters & p)
            : p(p)
        {
        }
        int clp2(int i) const
        {
            int e; std::frexp(i, &e);
            return 1 << e;
        }
        int sgn(double x) const
        {
            return (x < 0) ? -1 : (x > 0) ? 1 : 0;
        }
        void project(const bitmap & src, bitmap & dst)
        {
            w = src.w;
            h = src.h;
            S = std::hypot(src.w, src.h);
            C = clp2(S);

            dst.reshape(C, p.n);

            geom::point < int > p1, p2;

            double imax = 0;

            for (size_t t = 0; t < p.n; ++t)
            for (size_t i = 0; i < S; ++i)
            {
                int s = (int) i - (int) S / 2;
                make_line(w, h, s, M_PI / p.n * t, p1, p2);
                math::bresenham_rasterize(p1, p2, [&] (const geom::point < int > & p) {
                    dst.data[i][t] += src.data[p.y][p.x];
                });
                imax = max(imax, dst.data[i][t]);
            }

            for (size_t t = 0; t < p.n; ++t)
            for (size_t i = 0; i < S; ++i)
            {
                dst.data[i][t] /= imax;
            }
        }
        void backproject(const bitmap & src, bitmap & dst)
        {
            Q.clear();
            Q.resize(C);
            
            // calculate Q = /\/\ (regularized \/)
            for (size_t i = 0; i < S; ++i)
            {
                double s = ((int) i - (int) S / 2);
                double maxw = S / M_PI / 4;
                Q[i] = (sgn(2 * M_PI * maxw - s) + 2 * sgn(s) - sgn(2 * M_PI * maxw + s)) * std::sin(s / 2. / maxw);
            }

            bitmap prj;
            prj.reshape(S, p.n);

            // calculate filtered projection using FFT:
            // f = p ** q = FFT{ P * Q }
            for (size_t t = 0; t < p.n; ++t)
            {
                P.clear();
                P.resize(C);
                for (size_t i = 0; i < S; ++i)
                {
                    P[i] = src.data[i][t];
                }
                fourier(P.data(), C, -1);
                for (size_t i = 0; i < C; ++i)
                {
                    P[i] = P[i] * Q[i];
                }
                fourier(P.data(), C, 1);
                for (size_t i = 0; i < S; ++i)
                {
                    prj.data[i][t] = P[i].re;
                }
            }

            dst.reshape(h, w);

            double imax = 0;

            // calculate back projection using
            // simplest integration method
            // + without any interpolation of P[s,theta]
            for (int i = 0; i < h; ++i)
            for (int j = 0; j < w; ++j)
            for (size_t t = 0; t < p.n; ++t)
            {
                int s = (i - (int) h / 2) * std::cos(M_PI / p.n * t) - (j - (int) w / 2) * std::sin(M_PI / p.n * t) + S / 2;
                if (s < 0 || s >= S) continue;
                dst.data[i][j] += prj.data[s][t];
                imax = max(imax, dst.data[i][j]);
            }

            for (size_t i = 0; i < h; ++i)
            for (size_t j = 0; j < w; ++j)
            {
                dst.data[i][j] /= imax;
            }
        }
    };

    /*****************************************************/
    /*                     model                         */
    /*****************************************************/

    enum stage { stage_no, stage_source, stage_final };

    struct model_data
    {
        parameters params;
        model::stage stage;
        CBitmap csource;
        bitmap source;
        CBitmap cnoised;
        bitmap noised;
        CBitmap cprojection;
        bitmap projection;
        CBitmap crecovered;
        bitmap recovered;
    };
}