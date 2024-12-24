#pragma once
using namespace std;
#include "fixed.hpp"

#ifdef FLOAT
#error "FLOAT is already defined"
#endif
#ifdef DOUBLE
#error "DOUBLE is already defined"
#endif
#ifdef FIXED
#error "FIXED is already defined"
#endif
#ifdef FAST_FIXED
#error "FAST_FIXED is already defined"
#endif
#ifdef S
#error "S is defined"
#endif

#define FLOAT            float
#define DOUBLE           double
#define FAST_FIXED(N, K) types::FastFixed<N, K>
#define FIXED(N, K)      types::Fixed<N, K>

#define S(N, M)                  \
    types::SizePair {            \
        .rows = N, .columns = M, \
    }

namespace types {
struct data_from_file
{
    data_from_file(size_t N, size_t M, std::ifstream& in)
    {
        p.resize(N);
        velocity.resize(N);
        last_use.resize(N);
        for (int i = 0; i < N; i++)
        {
            p[i].resize(M);
            velocity[i].resize(M);
            last_use[i].resize(M);
        }
        in >> g;
        for (auto &elem : rho)
            in >> elem;
        in >> tick;
        for (auto &i : p)
            for (auto &j : i)
                in >> j;
        for (auto &i : velocity)
            for (auto &j : i)
                for (auto &k : j)
                    {in >> k;}
        for (auto &i : last_use)
            for (auto &j : i)
                {in >> j;}
        in >> UT;
    }
    double g;
    double rho[256];
    size_t tick = 0;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<std::array<double, 4>>> velocity;
    std::vector<std::vector<int>> last_use;
    int UT;  
};
struct SizePair {
    std::size_t rows{};
    std::size_t columns{};
};

struct Context {
    std::vector<std::vector<char>> field;
    data_from_file data;
};

template <class... Types>
struct TypesList;

namespace detail {

template <std::size_t I, class Type, class... Types>
struct TypesListGetAtHelper1 {
    static_assert(I < 1 + sizeof...(Types));
    using type = typename TypesListGetAtHelper1<I - 1, Types...>::type;
};

template <class Type, class... Types>
struct TypesListGetAtHelper1<0, Type, Types...> {
    using type = Type;
};

template <std::size_t I, class TypesListType>
struct TypesListGetAtHelper2;

template <std::size_t I, class... Types>
struct TypesListGetAtHelper2<I, TypesList<Types...>> {
    using type = typename TypesListGetAtHelper1<I, Types...>::type;
};

}  // namespace detail

template <std::size_t I, class TypesListType>
using TypesListGetAt = typename detail::TypesListGetAtHelper2<I, TypesListType>::type;

template <std::size_t I, class... Types>
using TypesPackGetAt = typename detail::TypesListGetAtHelper1<I, Types...>::type;

namespace detail {

inline constexpr std::size_t kDynamicExtent = std::numeric_limits<std::size_t>::max();

template <class PElementType,
          class VelocityElementType,
          class VelocityFlowElementType,
          std::size_t Rows    = kDynamicExtent,
          std::size_t Columns = kDynamicExtent>
class SimulationSession final {
    static_assert(Rows > 0);
    static_assert(Columns > 0);
    static_assert((Rows == kDynamicExtent && Columns == kDynamicExtent) ^
                  (Rows != kDynamicExtent && Columns != kDynamicExtent));

    static constexpr bool kUseStaticSize = Rows != kDynamicExtent && Columns != kDynamicExtent;

    static constexpr std::array<std::pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

    using RhoElementType = PElementType;

    template <class MatrixElementType>
    using StaticMatrixStorage = std::array<std::array<MatrixElementType, Columns>, Rows>;

    template <class MatrixElementType>
    using DynamicMatrixStorage = std::vector<std::vector<MatrixElementType>>;

    template <class MatrixElementType>
    using MatrixType = std::conditional_t<kUseStaticSize,
                                          StaticMatrixStorage<MatrixElementType>,
                                          DynamicMatrixStorage<MatrixElementType>>;

    template <class ElementType>
    struct VectorField {
        using Storage    = MatrixType<std::array<ElementType, deltas.size()>>;
        using value_type = Storage::value_type;

        Storage v{};
        void clear()
        {
            for (auto &i : v)
                for (auto &j : i)
                    for (auto &k : j)
                        k = ElementType(0);
        }
        VectorField()
            requires(kUseStaticSize)
        = default;

        VectorField(std::size_t dynamic_size, const value_type& row)
            requires(!kUseStaticSize)
            : v{dynamic_size, row} {}

        ElementType& add(int x, int y, int dx, int dy, const ElementType& dv) {
            return get(x, y, dx, dy) += dv;
        }

        ElementType& get(int x, int y, int dx, int dy) {
            size_t i = static_cast<std::size_t>(std::ranges::find(deltas, std::make_pair(dx, dy)) -
                                                deltas.begin());
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    using FieldCell           = char;
    using Field               = MatrixType<FieldCell>;
    using PStorage            = MatrixType<PElementType>;
    using VelocityStorage     = VectorField<VelocityElementType>;
    using VelocityFlowStorage = VectorField<VelocityFlowElementType>;
    using LastUseStorage      = MatrixType<int>;

public:
    explicit constexpr SimulationSession(const Context& ctx) noexcept
        requires(kUseStaticSize): saved(ctx.data)
    {
        assert(ctx.field.size() == Rows);
        assert(ctx.field.front().size() == Columns);

        for (std::size_t i = 0; auto& row : field) {
            const std::vector<char>& dynamic_field = ctx.field[i];
            assert(dynamic_field.size() == row.size());
            std::ranges::copy(dynamic_field, row.begin());
        }
    }

    explicit constexpr SimulationSession(const Context& ctx)
        requires(!kUseStaticSize)
        : SimulationSession(ctx.field, ctx.field.size(), ctx.field.front().size(), ctx.data) {}

    void start() {
        std::size_t rows    = field.size();
        std::size_t columns = field.front().size();
        assert(rows > 0);
        assert(columns > 0);
        std::cout << "Starting simulation with field size: (" << rows << ", " << columns <<")\n";
        int dirs[rows][columns]{};
        for (int i = 0; i < 256; i++)
        {
            rho[i] = (saved.rho[i]);
        }
        VelocityElementType g(saved.g);
        tick = saved.tick;
        if (tick != 0)
        {
            std::cout << "CHECK";
            UT = saved.UT;
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    p[i][j] = (saved.p[i][j]);
                }
            }
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    last_use[i][j] = (saved.last_use[i][j]);
                }
            }
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    velocity.v[i][j][0] = (saved.velocity[i][j][0]);
                    velocity.v[i][j][1] = (saved.velocity[i][j][1]);
                    velocity.v[i][j][2] = (saved.velocity[i][j][2]);
                    velocity.v[i][j][3] = (saved.velocity[i][j][3]);
                }
            }
        }
        for (size_t x = 0; x < rows; ++x) {
            for (size_t y = 0; y < columns; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy] : deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }
        for (size_t i = tick; i < 1000000; ++i) {
            PElementType total_delta_p(0);
            // Apply external forces
            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < columns; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
                }
            }

            // Apply forces from p
            for (int k = 0; k < rows; k++)
            {
                for (int j = 0; j < columns; j++)
                {
                    p_old[k][j] = p[k][j];
                }
            }
            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < columns; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && p_old[nx][ny] < p_old[x][y]) {
                            auto delta_p = p_old[x][y] - p_old[nx][ny];
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, VelocityElementType(force / rho[(int) field[x][y]]));
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow.clear();
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < columns; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < columns; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy] : deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = (VelocityElementType)new_v;
                            auto force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < rows; ++x) {
                for (size_t y = 0; y < columns; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        if (random01() < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                std::ofstream out("./local_save/save_"+std::to_string(rnd()));
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < columns; y++)
                    {
                        cout << field[x][y];
                        out << field[x][y];
                    }
                    out << endl;
                    cout << endl;
                }
                out << g << std::endl;
                for (int j = 0; j < 256; j++)
                {
                    out << rho[j] << " ";
                }
                out << endl;
                out << i << endl;
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < columns; y++)
                    {
                        out << p[x][y] << " ";
                    }
                    out << endl;
                }
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < columns; y++)
                    {
                        out << velocity.v[x][y][0] << " " << velocity.v[x][y][1] << " " << velocity.v[x][y][2] << " " << velocity.v[x][y][3] << " ";
                    }
                    out << endl;
                }
                for (size_t x = 0; x < rows; ++x) {
                    for (size_t y = 0; y < columns; y++)
                    {
                        out << last_use[x][y] << " ";
                    }
                    out << endl;
                }
                out << UT;
                out.close();
            }
        }
        assert(field.size() == rows);
        assert(p.size() == rows);
        assert(p_old.size() == rows);
        assert(velocity.v.size() == rows);
        assert(velocity_flow.v.size() == rows);
        assert(last_use.size() == rows);

        assert(field.front().size() == columns);
        assert(p.front().size() == columns);
        assert(p_old.front().size() == columns);
        assert(velocity.v.front().size() == columns);
        assert(velocity_flow.v.front().size() == columns);
        assert(last_use.front().size() == columns);

        assert(field.back().size() == columns);
        assert(p.back().size() == columns);
        assert(p_old.back().size() == columns);
        assert(velocity.v.back().size() == columns);
        assert(velocity_flow.v.back().size() == columns);
        assert(last_use.back().size() == columns);
    }
    tuple<PElementType, bool, std::pair<int, int>> propagate_flow(int x, int y, PElementType lim) { //CHECK
        last_use[x][y] = UT - 1;
        PElementType ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                PElementType rr = (PElementType)cap - (PElementType)flow;
                auto vp = min(lim, rr);
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, VelocityFlowElementType(vp));
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, (VelocityFlowElementType)t);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }
    PElementType random01() {
        return PElementType(double(rnd()) / rnd.max());
    }

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }
    PElementType move_prob(int x, int y) {
        PElementType sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    
    struct ParticleParams {
        
        char type;
        PElementType cur_p;
        array<VelocityElementType, deltas.size()> v;

        void swap_with(char &fieldElement, PElementType &p, array<VelocityElementType, deltas.size()> &velocity_element) {
            swap(fieldElement, type);
            swap(p, cur_p);
            swap(velocity_element, v);
        }
    };
    bool propagate_move(int x, int y, bool is_first) 
    {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            std::array<PElementType, deltas.size()> tres;
            PElementType sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            PElementType p = sum*random01();
            size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(field[x][y], p[x][y], velocity.v[x][y]);
                pp.swap_with(field[nx][ny], p[nx][ny], velocity.v[nx][ny]);
                pp.swap_with(field[x][y], p[x][y], velocity.v[x][y]);
            }
        }
        return ret;
    }


private:
    explicit constexpr SimulationSession(const Field& field, std::size_t rows, std::size_t columns, const data_from_file &ctx_data)
        requires(!kUseStaticSize)
        : field{field},
          p{rows, typename PStorage::value_type(columns)},
          p_old{rows, typename PStorage::value_type(columns)},
          velocity{rows, typename VelocityStorage::value_type(columns)},
          velocity_flow{rows, typename VelocityFlowStorage::value_type(columns)},
          last_use{rows, typename LastUseStorage::value_type(columns)},
          saved(ctx_data){}
    data_from_file saved;
    Field field{}; //save
    std::array<RhoElementType, 256> rho{}; //save
    PStorage p{}; //save
    PStorage p_old{};
    VelocityStorage velocity{}; //save
    VelocityFlowStorage velocity_flow{};
    LastUseStorage last_use{}; //save
    int UT = 0; //save
    size_t tick = 0; //save
    mt19937 rnd{1337};
};

template <class PElementType,
          class VelocityElementType,
          class VelocityFlowElementType,
          std::size_t Rows    = kDynamicExtent,
          std::size_t Columns = kDynamicExtent>
void start_simulation(const Context& ctx) {
    using SimulationSessionType = SimulationSession<PElementType, VelocityElementType,
                                                    VelocityFlowElementType, Rows, Columns>;
    SimulationSessionType{ctx}.start();
}

template <class PElementType,
          class VelocityElementType,
          class VelocityFlowElementType,
          SizePair StaticSize,
          SizePair... StaticSizes>
void select_size_and_start_simulation_impl(const Context& ctx) {
    static_assert(StaticSize.rows > 0);
    static_assert(StaticSize.columns > 0);

    if (StaticSize.rows == ctx.field.size() && StaticSize.columns == ctx.field.front().size()) {
        start_simulation<PElementType, VelocityElementType, VelocityFlowElementType,
                         StaticSize.rows, StaticSize.columns>(ctx);
    } else if constexpr (sizeof...(StaticSizes) == 0) {
        start_simulation<PElementType, VelocityElementType, VelocityFlowElementType>(ctx);
    } else {
        select_size_and_start_simulation_impl<PElementType, VelocityElementType,
                                              VelocityFlowElementType, StaticSizes...>(ctx);
    }
}

template <class PElementType,
          class VelocityElementType,
          class VelocityFlowElementType,
          SizePair... StaticSizes>
void select_size_and_start_simulation(const Context& ctx) {
    if constexpr (sizeof...(StaticSizes) > 0) {
        select_size_and_start_simulation_impl<PElementType, VelocityElementType,
                                              VelocityFlowElementType, StaticSizes...>(ctx);
    } else {
        start_simulation<PElementType, VelocityElementType, VelocityFlowElementType>(ctx);
    }
}

#define STRINGIFY_EXACT_NO_EVAL(expr) #expr

inline constexpr std::string_view kFloatTypeName   = STRINGIFY_EXACT_NO_EVAL(FLOAT);
inline constexpr std::string_view kDoubleTypeName  = STRINGIFY_EXACT_NO_EVAL(DOUBLE);
inline constexpr std::string_view kFastFixedPrefix = "FAST_FIXED(";
inline constexpr std::string_view kFastFixedSuffix = ")";
inline constexpr std::string_view kFixedPrefix     = "FIXED(";
inline constexpr std::string_view kFixedSuffix     = ")";

#undef STRINGIFY_EXACT_NO_EVAL

template <class AllowedFloatTypesList, class SelectedFloatTypesList, SizePair... StaticSizes>
class FloatTypesSelector;

template <class... AllowedFloatTypes, class... SelectedFloatTypes, SizePair... StaticSizes>
class FloatTypesSelector<TypesList<AllowedFloatTypes...>,
                         TypesList<SelectedFloatTypes...>,
                         StaticSizes...> {
    using FieldVec = std::vector<std::vector<char>>;

    static constexpr bool kCanUseFloatType =
        std::disjunction_v<std::is_same<AllowedFloatTypes, FLOAT>...>;

    static constexpr bool kCanUseDoubleType =
        std::disjunction_v<std::is_same<AllowedFloatTypes, DOUBLE>...>;

public:
    template <class... Args>
        requires std::conjunction_v<std::is_same<std::string_view, Args>...>
    static void select_type_and_size_and_start_simulation(const Context& ctx,
                                                          std::string_view type_name,
                                                          Args... type_names) {
        if (kCanUseFloatType && type_name == kFloatTypeName) {
            next_call_with_selected_type<FLOAT, Args...>(ctx, type_names...);
            return;
        }
        if (kCanUseDoubleType && type_name == kDoubleTypeName) {
            next_call_with_selected_type<DOUBLE, Args...>(ctx, type_names...);
            return;
        }

        std::size_t n{};
        std::size_t k{};
        if (try_parse_fixed_params(type_name, kFastFixedPrefix, kFastFixedSuffix, n, k)) {
            if (select_fixed_and_size_and_start_simulation<true, Args...>(n, k, ctx,
                                                                          type_names...)) {
                return;
            }
        }
        if (try_parse_fixed_params(type_name, kFixedPrefix, kFixedSuffix, n, k)) {
            if (select_fixed_and_size_and_start_simulation<false, Args...>(n, k, ctx,
                                                                           type_names...)) {
                return;
            }
        }

        throw_on_bad_type_name(type_name);
    }

private:
    static bool try_parse_fixed_params(std::string_view type_name,
                                       std::string_view prefix,
                                       std::string_view suffix,
                                       std::size_t& n,
                                       std::size_t& k) {
        if (!type_name.starts_with(prefix)) {
            return false;
        }
        type_name.remove_prefix(prefix.size());

        if (!type_name.ends_with(suffix)) {
            return false;
        }
        type_name.remove_suffix(suffix.size());

        const std::size_t sep_char_pos = type_name.find(',');
        if (sep_char_pos >= type_name.size()) {
            return false;
        }

        auto strip_sv = [](std::string_view s) noexcept {
            while (!s.empty() && std::isspace(s.front())) {
                s.remove_prefix(1);
            }
            while (!s.empty() && std::isspace(s.back())) {
                s.remove_suffix(1);
            }
            return s;
        };
        const std::string_view n_str = strip_sv(type_name.substr(0, sep_char_pos));
        const std::string_view k_str = strip_sv(type_name.substr(sep_char_pos + 1));

        if (std::from_chars(n_str.data(), n_str.data() + n_str.size(), n).ec != std::errc{}) {
            return false;
        }

        if (std::from_chars(k_str.data(), k_str.data() + k_str.size(), k).ec != std::errc{}) {
            return false;
        }

        return n > 0 && k > 0;
    }

    template <bool Fast, class... Args>
    static bool select_fixed_and_size_and_start_simulation(std::size_t n,
                                                           std::size_t k,
                                                           const Context& ctx,
                                                           Args... type_names) {
        return select_fixed_and_size_and_start_simulation_impl<0, Fast, Args...>(n, k, ctx,
                                                                                 type_names...);
    }

    template <std::size_t I, bool Fast, class... Args>
    static bool select_fixed_and_size_and_start_simulation_impl(std::size_t n,
                                                                std::size_t k,
                                                                const Context& ctx,
                                                                Args... type_names) {
        using FloatType = TypesPackGetAt<I, AllowedFloatTypes...>;

        if constexpr (requires {
                          {
                              FloatType::kNValue == std::size_t {}
                          } -> std::same_as<bool>;
                          {
                              FloatType::kKValue == std::size_t {}
                          } -> std::same_as<bool>;
                          {
                              FloatType::kFast == bool {}
                          } -> std::same_as<bool>;
                      }) {
            if constexpr (FloatType::kFast == Fast) {
                if (FloatType::kNValue == n && FloatType::kKValue == k) {
                    next_call_with_selected_type<FloatType, Args...>(ctx, type_names...);
                    return true;
                }
            }
        }

        if constexpr (I + 1 < sizeof...(AllowedFloatTypes)) {
            return select_fixed_and_size_and_start_simulation_impl<I + 1, Fast, Args...>(
                n, k, ctx, type_names...);
        }

        return false;
    }

    [[noreturn]] static void throw_on_bad_type_name(std::string_view type_name) {
        throw std::invalid_argument("Could not parse type " + std::string{type_name});
    }

    template <class FloatType, class... Args>
    static void next_call_with_selected_type(const Context& ctx, const Args&... type_names) {
        if constexpr (sizeof...(type_names) > 0) {
            using AllowedFloatTypesList  = TypesList<AllowedFloatTypes...>;
            using SelectedFloatTypesList = TypesList<SelectedFloatTypes..., FloatType>;

            FloatTypesSelector<AllowedFloatTypesList, SelectedFloatTypesList, StaticSizes...>::
                select_type_and_size_and_start_simulation(ctx, type_names...);
        } else {
            select_size_and_start_simulation<SelectedFloatTypes..., FloatType, StaticSizes...>(ctx);
        }
    }
};

}  // namespace detail

struct SimulationParams {
    std::string_view p_type_name{};
    std::string_view v_type_name{};
    std::string_view v_flow_type_name{};
};

template <class TypesList, SizePair... StaticSizes>
class Simulator;

template <class... FloatTypes, SizePair... StaticSizes>
class [[nodiscard]] Simulator<TypesList<FloatTypes...>, StaticSizes...> {
public:
    static Simulator from_params(SimulationParams params) {
        return Simulator{std::move(params)};
    }

    void start_on_field(const Context& ctx) const {
        assert(!ctx.field.empty());
        assert(!ctx.field.front().empty());
        detail::FloatTypesSelector<TypesList<FloatTypes...>, TypesList<>, StaticSizes...>::
            select_type_and_size_and_start_simulation(ctx, params_.p_type_name, params_.v_type_name,
                                                      params_.v_flow_type_name);
    }

private:
    explicit constexpr Simulator(SimulationParams&& params) noexcept : params_{std::move(params)} {}

    SimulationParams params_;
};

}  // namespace types
