#pragma once

namespace lazybastard {
/**
 * Class representing a toggle which can reach two possible states.
 *
 * A toggle has the size of a bool, but supports operations like an unsigned integer.
 *
 * Example:
 *
 * Operation: Multiplication
 *
 *    O1    O2    R
 * ___________________
 * |  1  |  1  |  1  |
 * |  1  | -1  | -1  |
 * | -1  | -1  |  1  |
 * | -1  |  1  | -1  |
 * -------------------
 *
 * All Toggle related operations:
 *
 *  Toggle(true)    &= true     ->  true
 *  Toggle(true)    &= false    ->  false
 *  Toggle(false)   &= false    ->  true
 *  Toggle(false)   &= true     ->  false
 *
 *  Toggle(true)    &= Toggle(true)     ->  true
 *  Toggle(true)    &= Toggle(false)    ->  false
 *  Toggle(false)   &= Toggle(false)    ->  true
 *  Toggle(false)   &= Toggle(true)     ->  false
 *
 *  Creation:
 *
 *  Toggle t = true;
 *
 *  Implicit conversion:
 *
 *  auto t = Toggle(false);
 *  t &= false; -> becomes true
 *
 *  if (t) {
 *      // t became true
 *  }
 */
struct Toggle {
  Toggle(bool b) : m_state(b){};
  operator bool() const { return m_state; };
  Toggle &operator&=(bool b) {
    if (m_state) {
      m_state &= b;
    } else {
      m_state = !m_state;
      m_state &= !b;
    }

    return *this;
  };
  Toggle operator&&(Toggle t) const { return m_state && t.m_state; };
  Toggle operator&&(bool b) const { return m_state && b; };
  Toggle operator==(Toggle t) const { return m_state == t.m_state; };
  Toggle operator!=(Toggle t) const { return !(*this == t); };

private:
  bool m_state{false};
};
} // namespace lazybastard