-- CNF Runner for NitroSat
local NitroSat = require("nitrosat")

-- Parse DIMACS CNF file
function parse_cnf(filename)
    local f = io.open(filename, "r")
    if not f then
        error("Cannot open file: " .. filename)
    end

    local num_vars = 0
    local num_clauses = 0
    local clauses = {}

    for line in f:lines() do
        -- Skip comments
        if line:sub(1, 1) == "c" then
            goto continue
        end

        -- Problem line
        if line:sub(1, 1) == "p" then
            local _, _, n, c = line:find("cnf%s+(%d+)%s+(%d+)")
            num_vars = tonumber(n)
            num_clauses = tonumber(c)
            goto continue
        end

        -- Clause line
        local literals = {}
        for lit in line:gmatch("(-?%d+)") do
            local l = tonumber(lit)
            if l ~= 0 then
                table.insert(literals, l)
            end
        end
        if #literals > 0 then
            table.insert(clauses, literals)
        end

        ::continue::
    end

    f:close()

    if num_vars == 0 then
        error("Invalid CNF file: no problem line found")
    end

    return {
        num_vars = num_vars,
        clauses = clauses
    }
end

-- Get filename from command line
local filename = arg[1]
if not filename then
    print("Usage: luajit run_cnf.lua <cnf_file>")
    print("Running on test.cnf as default...")
    filename = "/tmp/test.cnf"
end

print("Loading CNF from: " .. filename)
local instance = parse_cnf(filename)
print(string.format("Instance: %d variables, %d clauses", instance.num_vars, #instance.clauses))

-- Run solver
local solver = NitroSat.NitroSat.new(instance, { seed = 42 })
print("Running solver...")
local success, steps, sat_count = solver:solve({
    max_steps = 5000,
    lr = 0.002,
    verbose = true
})

print(string.format("Result: success=%s, steps=%d, satisfied=%d/%d",
    tostring(success), steps, sat_count, solver.num_clauses))
