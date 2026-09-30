"""Engineering requirements and explicit, immutable baseline comparisons."""

from __future__ import annotations

import math
import operator

OPERATORS = {"le": operator.le, "ge": operator.ge, "eq": operator.eq}


def evaluate(metrics: dict, requirements: list[dict]) -> dict:
    checks = []
    for rule in requirements:
        try:
            actual = metrics[rule["metric"]]["value"]
            if not math.isfinite(actual):
                raise ValueError("Nonfinite metric")
            status = "PASS" if OPERATORS[rule["operator"]](actual, rule["limit"]) else "FAIL"
            checks.append({**rule, "actual": actual, "status": status})
        except (KeyError, TypeError, ValueError) as exc:
            checks.append({**rule, "status": "ERROR", "reason": str(exc)})
    overall = (
        "NOT_EVALUATED"
        if not checks
        else "ERROR"
        if any(c["status"] == "ERROR" for c in checks)
        else "FAIL"
        if any(c["status"] == "FAIL" for c in checks)
        else "PASS"
    )
    return {"overall": overall, "checks": checks}


def regression(baseline: dict, candidate: dict, rules: list[dict]) -> dict:
    checks = []
    for rule in rules:
        try:
            if set(rule) != {"metric", "direction", "absolute", "relative"}:
                raise ValueError("Regression rule needs metric, direction, absolute, relative")
            if rule["direction"] not in ("lower", "higher"):
                raise ValueError("direction must be lower or higher")
            a, r = rule["absolute"], rule["relative"]
            if not all(
                isinstance(v, (int, float))
                and not isinstance(v, bool)
                and math.isfinite(v)
                and v >= 0
                for v in (a, r)
            ):
                raise ValueError("Allowances must be finite nonnegative numbers")
            b, c = baseline[rule["metric"]], candidate[rule["metric"]]
            if (b["unit"], b.get("version")) != (c["unit"], c.get("version")):
                raise ValueError("Metric unit/version mismatch")
            bv, cv = b["value"], c["value"]
            if not math.isfinite(bv) or not math.isfinite(cv):
                raise ValueError("Nonfinite metric")
            degradation = cv - bv if rule["direction"] == "lower" else bv - cv
            checks.append(
                {
                    **rule,
                    "baseline": bv,
                    "candidate": cv,
                    "allowance": a + r * abs(bv),
                    "status": "PASS" if degradation <= a + r * abs(bv) else "FAIL",
                }
            )
        except (KeyError, ValueError, TypeError) as exc:
            checks.append({**rule, "status": "ERROR", "reason": str(exc)})
    return {
        "overall": "ERROR"
        if not checks or any(x["status"] == "ERROR" for x in checks)
        else "FAIL"
        if any(x["status"] == "FAIL" for x in checks)
        else "PASS",
        "checks": checks,
    }
